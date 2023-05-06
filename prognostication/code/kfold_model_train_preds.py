import argparse
import warnings
import pandas as pd
import eli5
from eli5.sklearn import PermutationImportance
from eli5.permutation_importance import get_score_importances
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
from sksurv.metrics import integrated_brier_score, concordance_index_ipcw
from sksurv.linear_model import CoxnetSurvivalAnalysis
import torch
import torchtuples as tt
from pycox.evaluation import EvalSurv
from pycox.models import CoxTime
from pycox.models.cox_time import MLPVanillaCoxTime
from sklearn_pandas import DataFrameMapper
from sklearn.model_selection import RepeatedStratifiedKFold
import shap
from modules_preprocess import imputer_knn, data_normalize
from modules_mstate import *

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Cross validation on the validation cohort; helps to choose the model')
parser.add_argument('--seed', type=int, default=456, help='random seed (default: 42)')
parser.add_argument('--kfolds', type=int, default=3, help='Cross validation, number of folds')
parser.add_argument('--surv_model', type=str, default='neural_cox_non_prop', help='rsf, cph, neural_cox_non_prop, deep_surv')

parser.add_argument('--c_index_type', type=str, default='harrell', help='highly censored? | uno else harrell')
parser.add_argument('--n_repeat_kfold', type=int, default=1, help='times to repeat kfold')
parser.add_argument('--in_data_file', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in', help='input dataset path')
parser.add_argument('--path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1', help='input dataset path')
parser.add_argument('--legend', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/legend_PMMM_2023.xlsx', help='input dataset path')
parser.add_argument('--thresh_m_to_p2', type=int, default=365, help='days')

parser.add_argument('--perm_feat_imp', action='store_true', default=False, help='permutation feat imp')
parser.add_argument('--shap_feat_imp', action='store_true', default=False, help='SHAP features')
parser.add_argument('--get_feats_group', type=str, default='None', help='from master key file; if set None (100) then all feats')
parser.add_argument('--top_gen_feats', action='store_true', default=False, help='top features for genomics except rec. transloc.')
parser.add_argument('--top_n_gen_feats', type=int, default=1, help='top n genomic feats by perm importance')
parser.add_argument('--agg_reranked', action='store_true', default=False, help='based on - univariate & perm imp')

"""
Data prep
"""


def data_filter(df):
    df = df[df['duration'] >= 0]
    df = df.dropna(subset=['os_time', 'pfs_time', 'pfs_code', 'os_code'])

    df['time_SCT'].fillna(0, inplace=True)

    df_r_iss = pd.read_csv(args.in_data_file + '/r.iss.1933pts.txt', sep='\t')[['sample','R_ISS', 'Del_17p13.1']]
    df_r2_iss = pd.read_csv(args.in_data_file + '/r2.iss.1933pts.txt', sep='\t')[['sample', 'R2_ISS']]
    df_r_r2_iss = pd.merge(df_r_iss, df_r2_iss, on='sample', how='inner')
    df = pd.merge(df, df_r_r2_iss, on='sample', how='inner')

    # GEP70
    if args.get_feats_group == '10':
        df_gep = pd.read_csv(args.in_data_file + '/GEP70.txt', sep='\t')
        df = pd.merge(df, df_gep[['Patient', 'RNA_GEP70_dichotomous']], left_on='sample', right_on='Patient', how='inner')
        df.drop('Patient', axis=1, inplace=True)

    # map clinical vars
    map_gender = {'male': 0, 'female': 1}
    x_df = df.replace({'gender': map_gender})
    map_study = {'MMRF': 0, 'MPG': 1, 'Moffit': 2, 'MSKCC_292': 3, 'UAMS': 4}
    x_df = x_df.replace({'study': map_study})
    map_ecog = {'ecog<2': 0, 'ecog>=2': 1}
    x_df = x_df.replace({'ecog': map_ecog})
    map_iss = {'ISS1': 0, 'ISS2': 1, 'ISS3': 2}
    x_df = x_df.replace({'ISS': map_iss})
    map_r_iss = {'R-ISS1': 0, 'R-ISS2': 1, 'R-ISS3': 2}
    x_df = x_df.replace({'R_ISS': map_r_iss})
    map_r2_iss = {'I': 0, 'II': 1, 'III': 2, 'IV':3}
    x_df = x_df.replace({'R2_ISS': map_r2_iss})
    x_df["combo"] = np.select([x_df.combo <= 1, x_df.combo > 1], [0, 1])
    x_df["age"] = np.select([(x_df.age < 65), (x_df.age >= 65) & (x_df.age < 75), x_df.age >= 75], [0, 1, 2])

    map_ldh = {'Low': 0, 'High': 1, 'Normal': 2}
    x_df = x_df.replace({'LDH_level': map_ldh})

    map_race = {'WHITE': 0, 'BLACK OR AFRICAN AMERICAN': 1, 'OTHER': 2}
    x_df = x_df.replace({'Race': map_race})

    x_df['pfs_time'] = x_df['pfs_time'].astype(int)
    x_df['os_time'] = x_df['os_time'].astype(int)
    x_df['duration'] = x_df['duration'].astype(int)
    x_df['time_SCT'] = x_df['time_SCT'].astype(int)

    x_df['pfs_time'] = np.where(
        (x_df['duration'] > x_df['pfs_time']) & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
        x_df[['duration', 'pfs_time', 'os_time']].max().max(), x_df['pfs_time'])

    x_df['os_time'] = np.where(
        (x_df['duration'] > x_df['os_time']) & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
        x_df[['duration', 'pfs_time', 'os_time']].max().max(), x_df['os_time'])

    x_df['os_time'] = np.where(
        (x_df['pfs_time'] > x_df['os_time']) & (x_df['pfs_code'] == 1),
        x_df[['pfs_time', 'os_time']].max().max(), x_df['os_time'])

    threshold = args.thresh_m_to_p2
    x_df['duration'] = np.where(x_df['duration'] >= threshold, threshold, x_df['duration'])

    y_df = x_df[['pfs_code', 'pfs_time', 'os_code', 'os_time']]

    x_df['phases'] = np.where((x_df['phase'] == 'Phase1-only'), 'phase_1', 'phase_2')

    x_df.drop('phase', axis=1, inplace=True)

    x_df['stages'] = np.select([x_df['phases'] == 'phase_1',
                                   (x_df['phases'] == 'phase_2') & (x_df['pfs_time'] <= 540) & ~((x_df['pfs_code']==0) & (x_df['os_code']==1)),
                                   (x_df['phases'] == 'phase_2') & ((x_df['pfs_code']==0) & (x_df['os_code']==1)) & (x_df['os_time'] <= 540),
                                   (x_df['phases'] == 'phase_2') & (x_df['pfs_time'] > 540) & ~((x_df['pfs_code']==0) & (x_df['os_code']==1)),
                                   (x_df['phases'] == 'phase_2') & ((x_df['pfs_code']==0) & (x_df['os_code']==1)) & (x_df['os_time'] > 540)],
                                    ['induction', 'remission', 'remission', 'post-18months', 'post-18months'])

    x_df['sub_stages'] = np.select([(x_df['stages'] == 'induction') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'induction') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'induction') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'induction') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'remission') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1),
                                        (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 0),
                                        (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 0) & (x_df['os_code'] == 1),
                                       (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 0),
                                       (x_df['stages'] == 'post-18months') & (x_df['pfs_code'] == 1) & (x_df['os_code'] == 1)
                                       ],
                ['induction_00', 'induction_01', 'induction_10', 'induction_11',
                 'remission_00', 'remission_01', 'remission_10', 'remission_11',
                 'post-18months_00', 'post-18months_01', 'post-18months_10', 'post-18months_11'])

    df_x_entire = x_df.copy()

    x_df = x_df.drop(
        ['os_code', 'os_time', 'pfs_code', 'pfs_time', 'time_SCT', 'phases', 'stages', 'sub_stages', 'duration', 'study',
         'sample', 'DARA', 'ELO'], axis=1)

    # select groups of features for c-index
    x_df = select_features(x_df)
    # print('dataset shape is', x_df.shape)

    return x_df, y_df, df_x_entire


def select_features(df):
    df_fea_anno = pd.read_excel(args.legend)
    if args.get_feats_group == '10':
        df_fea_anno = df_fea_anno[
                            (df_fea_anno['column_id'] != 'time_SCT') & (df_fea_anno['column_id'] != 'study') &
                            (df_fea_anno['column_id'] != 'sample') & (df_fea_anno['column_id'] != 'phase') &
                            (df_fea_anno['column_id'] != 'os_time') & (df_fea_anno['column_id'] != 'pfs_time') &
                            (df_fea_anno['column_id'] != 'duration') & (df_fea_anno['column_id'] != 'pfs_code') &
                            (df_fea_anno['column_id'] != 'os_code') & (df_fea_anno['column_id'] != 'SCT_line')
                            & (df_fea_anno['column_id'] != 'DARA') & (df_fea_anno['column_id'] != 'ELO')
                            ]
    else:
        df_fea_anno = df_fea_anno[
            (df_fea_anno['column_id'] != 'time_SCT') & (df_fea_anno['column_id'] != 'study') &
            (df_fea_anno['column_id'] != 'sample') & (df_fea_anno['column_id'] != 'phase') &
            (df_fea_anno['column_id'] != 'os_time') & (df_fea_anno['column_id'] != 'pfs_time') &
            (df_fea_anno['column_id'] != 'duration') & (df_fea_anno['column_id'] != 'pfs_code') &
            (df_fea_anno['column_id'] != 'os_code') & (df_fea_anno['column_id'] != 'SCT_line')
            & (df_fea_anno['column_id'] != 'DARA') & (df_fea_anno['column_id'] != 'ELO') & (df_fea_anno['column_id'] != 'RNA_GEP70_dichotomous')
            ]

    if (args.get_feats_group == '23') or (args.get_feats_group == '24'): # r-iss / r2-iss split feature Del_17p13.1
        df_ = df
    else:
        df_ = df[df_fea_anno['column_id'].tolist()]

    if args.get_feats_group != 'None':
        master_df = pd.read_csv(args.path + '/feat_matrix/main_keeper.csv')
        df_ = df_[master_df.loc[int(args.get_feats_group), 'feat_combo'].split(' ')]

    if args.perm_feat_imp or args.shap_feat_imp or args.top_gen_feats or args.get_feats_group == 'None':
        df_ = df_.drop(['R_ISS', 'R2_ISS'], axis=1)

    return df_


def raw_data_read_prepare():
    df = pd.read_csv(args.in_data_file + '/PMMM_train_03302023.txt', sep='\t')

    df = df.drop('SCT_line', axis=1)

    df = df.sort_values(by=['sample'])

    x_df, y_df, df_x_entire = data_filter(df)

    return x_df, y_df, df_x_entire


"""
Models
"""


def cph(train_x, train_y, test_x, test_x_bystate, test_y_bystate):
    train_y = train_y[['binary_event', 'event_time']]
    test_y_bystate = test_y_bystate[['binary_event', 'event_time']]

    cph = CoxnetSurvivalAnalysis(l1_ratio=0.01, fit_baseline_model=True)

    train_x_std_norm, scaler = data_normalize(train_x, 'std_scaler')
    if scaler is not None:
        test_x_std_norm = scaler.transform(test_x)
        test_x_by_state_std_norm = scaler.transform(test_x_bystate)
    else:
        test_x_std_norm = test_x_bystate
        test_x_by_state_std_norm = test_x_std_norm

    ytrain = Surv.from_dataframe('binary_event', 'event_time', train_y)

    cph.fit(train_x_std_norm, ytrain)
    lower, upper = np.percentile(ytrain['event_time'], [10, 90])
    times_ytrain = np.arange(lower, upper)

    # train
    if args.c_index_type == 'harrell':
        c_ind_train = cph.score(train_x_std_norm, ytrain)
    elif args.c_index_type == 'uno':
        try:
            surv_prob_tr = [fn for fn in cph.predict(train_x_std_norm)]
            c_score_uno = concordance_index_ipcw(ytrain, ytrain, surv_prob_tr)
            c_ind_train = c_score_uno[0]
        except Exception:
            c_ind_train = None
            pass

    # val by state
    yval = Surv.from_dataframe('binary_event', 'event_time', test_y_bystate)
    if args.c_index_type == 'harrell':
        c_ind_val = cph.score(test_x_by_state_std_norm, yval)
    elif args.c_index_type == 'uno':
        try:
            surv_prob_val = [fn for fn in cph.predict(test_x_by_state_std_norm)]
            c_score_uno = concordance_index_ipcw(ytrain, yval, surv_prob_val)
            c_ind_val = c_score_uno[0]
        except Exception:
            c_ind_val = None
            pass

    surv_tr_prob = np.row_stack([fn(times_ytrain) for fn in cph.predict_survival_function(train_x_std_norm)])
    surv_train = pd.DataFrame(surv_tr_prob.T, index=times_ytrain).astype('float32')
    surv_train.index.names = ['duration']

    surv_val_prob = np.row_stack([fn(times_ytrain) for fn in cph.predict_survival_function(test_x_std_norm)])
    surv_val = pd.DataFrame(surv_val_prob.T, index=times_ytrain).astype('float32')
    surv_val.index.names = ['duration']

    return surv_val, c_ind_train, c_ind_val, cph


def rsf(train_x, train_y, test_x, test_x_bystate, test_y_bystate, model_param):
    train_y = train_y[['binary_event', 'event_time']]
    test_y_bystate = test_y_bystate[['binary_event', 'event_time']]

    model = RandomSurvivalForest(n_estimators=model_param, oob_score=True, random_state=args.seed, min_samples_split=2,
                               min_samples_leaf=1, max_features='sqrt', max_depth=30)

    model.fit(train_x, Surv.from_dataframe('binary_event', 'event_time', train_y))

    surv_train = model.predict_survival_function(train_x, return_array=True)
    surv_train = pd.DataFrame(surv_train.T)
    surv_train.index = model.event_times_
    surv_train.index.names = ['duration']

    surv_val = model.predict_survival_function(test_x, return_array=True)
    surv_val = pd.DataFrame(surv_val.T)
    surv_val.index = model.event_times_
    surv_val.index.names = ['duration']

    surv_val_bystate = model.predict_survival_function(test_x_bystate, return_array=True)
    surv_val_bystate = pd.DataFrame(surv_val_bystate.T)
    surv_val_bystate.index = model.event_times_
    surv_val_bystate.index.names = ['duration']

    ytrain = Surv.from_dataframe('binary_event', 'event_time', train_y)
    # train
    if args.c_index_type == 'harrell':
        c_ind_train = model.score(train_x, ytrain)
    elif args.c_index_type == 'uno':
        try:
            risk_score = model.predict(train_x)
            c_score_uno = concordance_index_ipcw(ytrain, ytrain, risk_score)
            c_ind_train = c_score_uno[0]
        except Exception:
            c_ind_train = None
            pass

    # val by state
    yval = Surv.from_dataframe('binary_event', 'event_time', test_y_bystate)
    if args.c_index_type == 'harrell':
        c_ind_val = model.score(test_x_bystate, yval)
    elif args.c_index_type == 'uno':
        try:
            risk_score = model.predict(test_x)
            c_score_uno = concordance_index_ipcw(ytrain, yval, risk_score)
            c_ind_val = c_score_uno[0]
        except Exception:
            c_ind_val = None
            pass

    return surv_val, c_ind_train, c_ind_val, model


def cox_time(train_x, train_y, test_x, test_x_bystate, test_y_bystate, cols_x, model_param, group, df_test):
    np.random.seed(args.seed)
    _ = torch.manual_seed(model_param[2])

    cols_leave = [x for x in cols_x ]
    leave = [(col, None) for col in cols_leave]
    x_mapper = DataFrameMapper(leave)

    train_x = pd.DataFrame(train_x, columns=cols_x)
    train_x = x_mapper.fit_transform(train_x).astype('float32')

    test_x = pd.DataFrame(test_x, columns=cols_x)
    test_x = x_mapper.transform(test_x).astype('float32')

    test_x_bystate = pd.DataFrame(test_x_bystate, columns=cols_x)
    test_x_bystate = x_mapper.transform(test_x_bystate).astype('float32')

    train_x = train_x[:, :len(cols_x)]
    test_x = test_x[:, :len(cols_x)]
    test_x_bystate = test_x_bystate[:, :len(cols_x)]

    train_y = train_y[['binary_event', 'event_time']]
    test_y_bystate = test_y_bystate[['binary_event', 'event_time']]

    labtrans = CoxTime.label_transform()
    get_target = lambda df: (df['event_time'].values, df['binary_event'].values)

    y_train = labtrans.fit_transform(*get_target(train_y))
    y_val = labtrans.transform(*get_target(test_y_bystate))
    val = tt.tuplefy(test_x_bystate, y_val)

    in_features = train_x.shape[1]
    num_nodes = model_param[0]
    batch_norm = False
    dropout = model_param[4]
    net = MLPVanillaCoxTime(in_features, num_nodes, batch_norm, dropout)

    optimizer = tt.optim.AdamWR(decoupled_weight_decay=0.025, cycle_eta_multiplier=0.5)
    model = CoxTime(net, optimizer, labtrans=labtrans)

    model.optimizer.set_lr(0.005)
    batch_size = model_param[1]
    epochs = 20
    callbacks = [tt.cb.EarlyStopping()]
    log = model.fit(train_x, y_train, batch_size, epochs, callbacks, val_data=val, val_batch_size=batch_size, verbose=False)
    _ = log.plot()

    _ = model.compute_baseline_hazards()
    surv_train = model.predict_surv_df(train_x)
    surv_val = model.predict_surv_df(test_x)
    surv_val_bystate = model.predict_surv_df(test_x_bystate)

    # train
    durations_tr, events_tr = get_target(train_y)
    ev = EvalSurv(surv_train, durations_tr, events_tr, censor_surv='km')
    c_ind_train = ev.concordance_td('antolini')

    # val
    durations_val, events_val = get_target(test_y_bystate)
    ev = EvalSurv(surv_val_bystate, durations_val, events_val, censor_surv='km')
    c_ind_val = ev.concordance_td('antolini')

    df_test_x_bystate = pd.DataFrame(test_x_bystate, columns=cols_x)

    # construct shap values
    if args.shap_feat_imp and group != 'Move to phase 2' and group != 'Non-progress & deceased (P2)':
        explainer = shap.KernelExplainer(model.predict_surv, shap.kmeans(train_x, 2))
        shap_values = explainer.shap_values(df_test_x_bystate)

        mean_shap_values = np.mean(shap_values, axis=0)
        df_shap_mean = pd.DataFrame(mean_shap_values, columns=cols_x)
        df_test_shap = df_test_x_bystate
        print('shap shape', df_test_shap.shape, df_test.shape)
        #df_test_shap.index = df_test['sample']
    else:
        df_shap_mean = None
        df_test_shap = None

    return surv_val, c_ind_train, c_ind_val, model, df_shap_mean, df_test_shap


def nnet_survival_score(X, y):
    surv = y[0].predict_surv_df(X.astype('float32'))
    get_target = lambda df: (df['event_time'], df['binary_event'])
    durations_val, events_val = get_target(y[1])
    ev = EvalSurv(surv, durations_val, events_val, censor_surv='km')
    c_ind = ev.concordance_td('antolini')
    return c_ind


"""
Main caller to train, predict and compute feature importances
"""


def train_survival_(fold, train_x, train_y, test_x, test_x_bystate, test_y_bystate, df_test, group, model_param):
    # print('train surv', group)
    overall_fea_class = pd.read_excel(args.legend)

    # aggregated column wise pred probs
    df_pred_probs = pd.DataFrame()

    # feat imps
    df_feat_imp = pd.DataFrame()
    cols_x = list(train_x)

    # impute
    train_x, imput_model = imputer_knn(train_x, 'train', '')
    test_x, _ = imputer_knn(test_x, 'val', imput_model)
    test_x_bystate, _ = imputer_knn(test_x_bystate, 'val', imput_model)

    df_train_y = train_y.copy();  df_test_y_bystate = test_y_bystate.copy()

    if args.surv_model == 'cph':
        val_surv_df, c_ind_train, c_ind_val, estimator = cph(train_x, df_train_y, test_x, test_x_bystate, df_test_y_bystate)

    if args.surv_model == 'rsf':
        val_surv_df, c_ind_train, c_ind_val, estimator = rsf(train_x, df_train_y, test_x,
                                                                    test_x_bystate, df_test_y_bystate, model_param)

    if args.surv_model == 'neural_cox_non_prop':
        val_surv_df, c_ind_train, c_ind_val, estimator, mean_shap_values, df_test_shap = \
                                                    cox_time(train_x, df_train_y, test_x,
                                                                test_x_bystate, df_test_y_bystate, cols_x, model_param,
                                                                                 group, df_test)

    # interpolate
    idx = np.arange(0, 2554, 1)
    df_interp_pred = pd.merge(val_surv_df, pd.DataFrame(index=idx), left_index=True, right_index=True, how='outer')
    df_interp_pred = df_interp_pred.interpolate()
    df_interp_pred = df_interp_pred.head(2554)

    df_pred_probs = pd.concat((df_pred_probs, df_interp_pred), axis=1)

    if group == 'Move to phase 2':
        df_pred_probs.iloc[args.thresh_m_to_p2+1:] = 0

    df_pred_probs = df_pred_probs.fillna(1)

    df_pred_probs.columns = df_test['sample'].tolist()

    ##### prep training set for permutation importance
    train_y = Surv.from_dataframe('binary_event', 'event_time', df_train_y[['binary_event', 'event_time']])

    # permutation importance rsf
    if args.perm_feat_imp == True and args.surv_model == 'rsf':
        perm = PermutationImportance(estimator, n_iter=5, random_state=args.seed)
        perm.fit(train_x, train_y)
        eli_feat_imp = eli5.format_as_dataframe(eli5.explain_weights(perm, top=None))
        eli_feat_imp['feature_ind'] = eli_feat_imp['feature'].str.replace(r'x', '')
        eli_feat_imp.drop('feature', axis=1, inplace=True)
        eli_feat_imp['feature_ind'] = eli_feat_imp['feature_ind'].astype(int)
        eli_feat_imp['feature'] = eli_feat_imp['feature_ind'].apply(lambda i: cols_x[i])
        eli_feat_imp = eli_feat_imp.sort_values('feature')
        df_feat_imp = pd.concat([df_feat_imp, eli_feat_imp], ignore_index=True)

    # permutation feature importance deep nets
    if args.perm_feat_imp == True and ('cph' not in args.surv_model or args.surv_model != 'rsf'):
        base_score, score_decreases = get_score_importances(nnet_survival_score, train_x, (estimator, train_y), n_iter=5)
        feature_importances = np.mean(score_decreases, axis=0)
        eli_df = pd.DataFrame(feature_importances, index=cols_x).reset_index().rename(
            columns={'index': 'feature', 0: 'weight'})
        eli_df['std'] = np.std(score_decreases, axis=0)
        df_feat_imp = pd.concat([df_feat_imp, eli_df], ignore_index=True)

    # merge with legend
    if args.perm_feat_imp == True and ('cph' not in args.surv_model):
        df_feat_imp = pd.merge(df_feat_imp, overall_fea_class[['column_id', 'display_id', 'class_generic', 'class_granular']],
                               left_on='feature', right_on='column_id', how='left')
        df_feat_imp['group'] = group

    if not (args.shap_feat_imp and args.surv_model == 'neural_cox_non_prop'):
        mean_shap_values = None
        df_test_shap = None

    df_metrics = pd.DataFrame({'C_index_train': [c_ind_train], 'C_index_val': [c_ind_val], 'Kfolds': [fold],
                               'Fold': [args.kfolds], 'Model': [args.surv_model], 'Group': [group]})

    return df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def get_top_features(t_x, gp):
    if args.agg_reranked:
        df_feat_rank = pd.read_csv(args.path + '/permut_feat_imps/aggregate_reranked/' + gp +
                                   '_' + args.surv_model + '_folds_feat_imp.csv')
        print('re-rank-feats in main file')
    else:
        df_feat_rank = pd.read_csv(args.path + '/permut_feat_imps/aggregate/' + gp +
                                    '_' + args.surv_model + '_folds_feat_imp.csv')

    overall_fea_class = pd.read_excel(args.legend)

    df_feat_rank = pd.merge(df_feat_rank, overall_fea_class[['column_id', 'class_generic', 'class_granular']], on='column_id', how='inner')

    # baseline feats = ISS + Clinical + Treatment + Cont. treat + SCT + Rec. Translocation
    baseline_feats = df_feat_rank[(df_feat_rank['class_generic'] == 'ISS') | (df_feat_rank['class_generic'] == 'Treatment')
                              | (df_feat_rank['class_generic'] == 'Clinical') | (df_feat_rank['class_generic'] == 'Continuous treatment')
                              | (df_feat_rank['class_generic'] == 'SCT')
                                  #| (df_feat_rank['class_granular'] == 'Translocation')
                            ]

    # to get top n genomic feats - CNV/Mutations
    top_gen_feats = df_feat_rank[(df_feat_rank['class_granular'] == 'CNV') | (df_feat_rank['class_granular'] == 'CNV_Mutations')
                                 | (df_feat_rank['class_granular'] == 'Mutations')].iloc[:args.top_n_gen_feats]

    # baseline (with transloc) + top genomic feats + top univariate feats (intersection of os and pfs)
    if args.agg_reranked:
        uni_imp_feats = pd.read_csv(args.path + '/plots/df_feat_uni.csv')
        cols_subset = uni_imp_feats['Genomic.Feature'].tolist() + top_gen_feats['column_id'].tolist() + baseline_feats['column_id'].tolist()
    else:
        cols_subset = top_gen_feats['column_id'].tolist() + baseline_feats['column_id'].tolist()

    t_x = t_x[cols_subset]
    return t_x


"""
Multistate methods
"""


def progress_p1(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list):
    group = 'Progress (P1)'

    tra_x_bystate, tra_y_bystate = crit_prog_p1(train_x, train_y, df_train)
    tes_x_bystate, tes_y_bystate = crit_prog_p1(test_x, test_y, df_test)

    if args.top_gen_feats:
        tra_x_bystate = get_top_features(tra_x_bystate, group)
        tes_x_bystate = get_top_features(tes_x_bystate, group)
        test_x = get_top_features(test_x, group)

    tra_x_bystate_cp, tes_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), tes_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        tes_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, test_x_cp, tes_x_bystate_cp, tes_y_bystate, df_test, group, model_param_list[0])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_dec_p1(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list):
    group = 'Progress & deceased (P1)'

    tra_x_bystate, tra_y_bystate = crit_prog_dec_p1(train_x, train_y, df_train)
    tes_x_bystate, tes_y_bystate = crit_prog_dec_p1(test_x, test_y, df_test)

    if args.top_gen_feats:
        tra_x_bystate = get_top_features(tra_x_bystate, group)
        tes_x_bystate = get_top_features(tes_x_bystate, group)
        test_x = get_top_features(test_x, group)

    tra_x_bystate_cp, tes_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), tes_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        tes_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, test_x_cp,
            tes_x_bystate_cp, tes_y_bystate, df_test, group, model_param_list[1])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def move_to_phase2(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list):
    group = 'Move to phase 2'

    tra_x_bystate, tra_y_bystate = crit_move_to_p2(train_x, train_y, df_train)
    tes_x_bystate, tes_y_bystate = crit_move_to_p2(test_x, test_y, df_test)

    if args.top_gen_feats:
        tra_x_bystate = get_top_features(tra_x_bystate, group)
        tes_x_bystate = get_top_features(tes_x_bystate, group)
        test_x = get_top_features(test_x, group)

    tra_x_bystate_cp, tes_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), tes_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        tes_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, test_x_cp,
            tes_x_bystate_cp, tes_y_bystate, df_test, group, model_param_list[2])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_p2(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list):
    group = 'Progress (P2)'

    tra_x_bystate, tra_y_bystate = crit_prog_p2(train_x, train_y, df_train)
    tes_x_bystate, tes_y_bystate = crit_prog_p2(test_x, test_y, df_test)

    if args.top_gen_feats:
        tra_x_bystate = get_top_features(tra_x_bystate, group)
        tes_x_bystate = get_top_features(tes_x_bystate, group)
        test_x = get_top_features(test_x, group)

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate, tra_y_bystate, test_x, tes_x_bystate, tes_y_bystate,
                                df_test, group, model_param_list[3])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_dec_p2(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list):
    group = 'Progress & deceased (P2)'

    tra_x_bystate, tra_y_bystate = crit_prog_dec_p2(train_x, train_y, df_train)
    tes_x_bystate, tes_y_bystate = crit_prog_dec_p2(test_x, test_y, df_test)

    if args.top_gen_feats:
        tra_x_bystate = get_top_features(tra_x_bystate, group)
        tes_x_bystate = get_top_features(tes_x_bystate, group)
        test_x = get_top_features(test_x, group)

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate, tra_y_bystate, test_x, tes_x_bystate, tes_y_bystate,
                                df_test, group, model_param_list[4])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def non_progress_dec_p2(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list):
    group = 'Non-progress & deceased (P2)'

    tra_x_bystate, tra_y_bystate = crit_non_prog_dec_p2(train_x, train_y, df_train)
    tes_x_bystate, tes_y_bystate = crit_non_prog_dec_p2(test_x, test_y, df_test)

    if args.top_gen_feats:
        tra_x_bystate = get_top_features(tra_x_bystate, group)
        tes_x_bystate = get_top_features(tes_x_bystate, group)
        test_x = get_top_features(test_x, group)

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate, tra_y_bystate, test_x, tes_x_bystate, tes_y_bystate,
                                df_test, group, model_param_list[5])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def phase_1(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap):
    # print('phase 11', train_y['os_time'].min(), train_y['os_time'].max())

    # G1 Phase 1 induction (progressors vs non-progressors)
    df_pred_probs_prog_p1, gp_prog_p1, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_p1(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list)
    dict_pred_probs[gp_prog_p1] = df_pred_probs_prog_p1
    dict_feats[gp_prog_p1] = df_feat_imp
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_mean_shap_values[gp_prog_p1] = mean_shap_values
    dict_concat_test_shap[gp_prog_p1] = df_test_shap
    # print('phase 12', train_y['os_time'].min(), train_y['os_time'].max())

    # G2 Phase 1 induction progressors alive vs deceased
    df_pred_probs_prog_dec_p1, gp_prog_dec_p1, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_dec_p1(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list)
    dict_pred_probs[gp_prog_dec_p1] = df_pred_probs_prog_dec_p1
    dict_feats[gp_prog_dec_p1] = df_feat_imp
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_mean_shap_values[gp_prog_dec_p1] = mean_shap_values
    dict_concat_test_shap[gp_prog_dec_p1] = df_test_shap
    # print('phase 13', train_y['os_time'].min(), train_y['os_time'].max())

    # G2 phase 2 vs 1
    df_pred_probs_to_p2, gp_to_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = move_to_phase2(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list)
    dict_pred_probs[gp_to_p2] = df_pred_probs_to_p2
    dict_feats[gp_to_p2] = df_feat_imp
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_mean_shap_values[gp_to_p2] = mean_shap_values
    dict_concat_test_shap[gp_to_p2] = df_test_shap

    return dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap


def phase_2(fold, train_x, train_y, test_x, test_y, df_train, df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap):
    df_p2 = df_train[(df_train['phases'] == 'phase_2')]
    X_p2 = train_x[train_x.index.isin(df_p2.index)]
    Y_p2 = train_y[train_y.index.isin(df_p2.index)]

    # reset
    post_induction = df_p2[(df_p2['phases'] == 'phase_2')]
    Y_p2.loc[post_induction.index.tolist(), 'os_time'] = Y_p2.loc[post_induction.index.tolist(), 'os_time'] - post_induction['duration']
    Y_p2.loc[post_induction.index.tolist(), 'pfs_time'] = Y_p2.loc[post_induction.index.tolist(), 'pfs_time'] - \
                                                             post_induction['duration']

    X_p2 = X_p2.reset_index(drop=True)
    Y_p2 = Y_p2.reset_index(drop=True)
    df_p2 = df_p2.reset_index(drop=True)
    #print('phase 21', Y_p2['os_time'].min(), Y_p2['os_time'].max())

    # G3 Phase 2 progress vs non-progress
    df_pred_probs_prog_p2, gp_prog_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_p2(fold, X_p2, Y_p2, test_x, test_y, df_p2, df_test, model_param_list)
    dict_pred_probs[gp_prog_p2] = df_pred_probs_prog_p2
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_feats[gp_prog_p2] = df_feat_imp
    dict_mean_shap_values[gp_prog_p2] = mean_shap_values
    dict_concat_test_shap[gp_prog_p2] = df_test_shap
    #print('phase 22', Y_p2['os_time'].min(), Y_p2['os_time'].max())

    # G5 Phase 2 progress alive vs progress dec
    df_pred_probs_prog_dec_p2, gp_prog_dec_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_dec_p2(fold, X_p2, Y_p2, test_x, test_y, df_p2, df_test, model_param_list)
    dict_pred_probs[gp_prog_dec_p2] = df_pred_probs_prog_dec_p2
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_feats[gp_prog_dec_p2] = df_feat_imp
    dict_mean_shap_values[gp_prog_dec_p2] = mean_shap_values
    dict_concat_test_shap[gp_prog_dec_p2] = df_test_shap
    # print('phase 23', Y_p2['os_time'].min(), Y_p2['os_time'].max())

    # G4 Phase 2 non-prog dec vs rest in p2
    df_pred_probs_non_prog_dec_p2, gp_non_prog_dec_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = non_progress_dec_p2(fold, X_p2, Y_p2, test_x, test_y, df_p2, df_test, model_param_list)
    dict_pred_probs[gp_non_prog_dec_p2] = df_pred_probs_non_prog_dec_p2
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_feats[gp_non_prog_dec_p2] = df_feat_imp
    dict_mean_shap_values[gp_non_prog_dec_p2] = mean_shap_values
    dict_concat_test_shap[gp_non_prog_dec_p2] = df_test_shap

    return dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap


"""
Stratified CV caller. This is the core method to call all the processes
"""
def patient_cross_val(model_param_list):
    X, Y, x_df = raw_data_read_prepare()

    # factorize to maintain proportion for stratified CV
    x_df_copy = x_df.copy()
    x_df_copy = x_df_copy.fillna(1000)
    x_df_copy['unique_patient_groups'] = pd.factorize(pd._libs.lib.fast_zip([x_df_copy['ISS'].values,
                                                x_df_copy['age'].values, x_df_copy['SCT_first_line'].values]))[0]

    df_concat_metrics = pd.DataFrame()
    dict_pred_probs = {}
    dict_feats = {}
    dict_mean_shap_values = {}
    dict_concat_test_shap = {}

    rskf = RepeatedStratifiedKFold(n_splits=args.kfolds, n_repeats=args.n_repeat_kfold, random_state=np.random.RandomState(2652124))
    fold = 0
    for train_index, test_index in rskf.split(X, x_df_copy['unique_patient_groups']):
        # print(fold, train_index.shape, test_index.shape)

        train_x, test_x = X.iloc[train_index].reset_index(drop=True), X.iloc[test_index].reset_index(drop=True)
        train_y, test_y = Y.iloc[train_index].reset_index(drop=True), Y.iloc[test_index].reset_index(drop=True)
        x_df_train = x_df.iloc[train_index].reset_index(drop=True)
        x_df_test = x_df.iloc[test_index].reset_index(drop=True)

        ########## Phase 1 groups ################
        dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap = phase_1(fold, train_x, train_y, test_x, test_y, x_df_train,
                            x_df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap)

        ########## Phase 2 groups ################
        dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap = phase_2(fold, train_x, train_y, test_x, test_y, x_df_train,
                            x_df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap)

        # model probs
        for ind, (k, v) in enumerate(dict_pred_probs.items()):
            if (args.get_feats_group != 'None' or args.get_feats_group == 'None') and args.top_gen_feats == False \
                                                                and not args.perm_feat_imp and not args.shap_feat_imp:
                if args.get_feats_group == 'None':
                    group_id = '100'
                else:
                    group_id = args.get_feats_group

                v.to_csv(args.path + '/kfold/model_probs/group_' + group_id + '_fold_' + str(fold) + '_' +
                    args.surv_model + '_' + k + '.csv', index=False)

        # shap values
        if args.shap_feat_imp:
            for ind, (k, v) in enumerate(dict_mean_shap_values.items()):
                if v is not None:
                    v.to_csv(args.path + '/shap_feat_imps/fold_' + str(fold) + '_' +
                                                k + '_' + args.surv_model + '_' + str(args.kfolds) + '.csv', index=False)

            for ind, (k, v) in enumerate(dict_concat_test_shap.items()):
                if v is not None:
                    v.to_csv(args.path + '/shap_feat_imps/test_data_fold_' + str(fold) + '_' +
                                                k + '_' + args.surv_model + '_' + str(args.kfolds) + '.csv', index=True)

        # permutation importance
        if args.perm_feat_imp:
            for ind, (k, feat_imp) in enumerate(dict_feats.items()):
                feat_imp[['column_id', 'display_id', 'class_generic', 'class_granular', 'weight', 'std', 'group']]\
                        .to_csv(args.path + '/permut_feat_imps/compose/fold_' + str(fold) + '_all_feat_weights_' +
                                k + '_' + args.surv_model + '_' + str(args.kfolds) + '.csv', index=False)

        fold += 1

    # metrics by state models
    if (args.get_feats_group != 'None' or args.get_feats_group == 'None') and args.top_gen_feats == False \
                                                                and not args.perm_feat_imp and not args.shap_feat_imp:
        # with subset of feats / with all feats
        if args.get_feats_group == 'None':
            group_id = '100'
        else:
            group_id = args.get_feats_group
        data_out_path = args.path + '/kfold/metrics/model_feat_combo/'
        feats_group = group_id

    elif args.top_gen_feats == True:
        data_out_path = args.path + '/kfold/metrics/model_genomics/'
        feats_group = str(args.top_n_gen_feats)

    if not args.perm_feat_imp and not args.shap_feat_imp:
        df_concat_metrics.to_csv(data_out_path + 'surv_group_' + feats_group + '_' +
                             args.surv_model + '_' + str(args.kfolds) + 'folds.csv',
                             index=False)

    return None


def main():
    """
    Script to build a model to predict progression or no progression
    :return: None
    """
    global args
    args = parser.parse_args()

    # print('algorithm is', args.surv_model)
    model_param = {}; model_param_list = []

    if args.surv_model == 'rsf':
        model_param_list = [750, 400, 750, 750, 750, 750]

    if args.surv_model == 'neural_cox_non_prop':
        model_param['num_nodes'] = [[128, 128], [128, 96], [64, 64], [128, 96], [128, 128], [98, 36]]
        model_param['batch_size'] = [128, 64, 128, 128, 128, 128]
        model_param['seed'] = [666, 444, 666, 787, 888, 888]
        model_param['rate_optim'] = [0.01, 0.001, 0.01, 0.01, 0.01, 0.01]
        model_param['dropout'] = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        model_param_list = [(i, j, k, l, m) for i, j, k, l, m in zip(model_param['num_nodes'],
                                        model_param['batch_size'], model_param['seed'],
                                        model_param['rate_optim'], model_param['dropout'] )]

    if 'cph' in args.surv_model:
        model_param_list = [None]*6

    print('start')

    patient_cross_val(model_param_list)

    print('done')


if __name__ == '__main__':
    main()