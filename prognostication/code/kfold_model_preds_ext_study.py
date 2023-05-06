# Prediction 1.0
# @ Arjun Raj Rajanna

import pandas as pd
import numpy as np
import seaborn as sn
import argparse
import warnings
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
from modules_preprocess import imputer_knn, data_normalize
from modules_mstate import *

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='CV model testing: inference of each model on external test cohort')
parser.add_argument('--seed', type=int, default=456, help='random seed (default: 42)')
parser.add_argument('--kfolds', type=int, default=3, help='Cross validation, number of folds')
parser.add_argument('--surv_model', type=str, default='neural_cox_non_prop', help='rsf, cph, neural_cox_non_prop, deep_surv')
parser.add_argument('--c_index_type', type=str, default='harrell', help='highly censored? | uno else harrell')
parser.add_argument('--n_repeat_kfold', type=int, default=1, help='times to repeat kfold')
parser.add_argument('--in_data_file', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in', help='input dataset path')
parser.add_argument('--path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1', help='input dataset path')
parser.add_argument('--legend', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/legend_PMMM_2023.xlsx', help='input dataset path')
parser.add_argument('--thresh_m_to_p2', type=int, default=365, help='days')
parser.add_argument('--get_feats_group', type=str, default='None', help='101 is test set ; from master key file; if None (100) then all feats')


def data_filter(df, ext_val=False):
    df = df[df['duration'] >= 0]
    df = df.dropna(subset=['os_time', 'pfs_time', 'pfs_code', 'os_code'])

    df['time_SCT'].fillna(0, inplace=True)

    if ext_val is not True: # not for independent Heidelberg cohort
        df_r_iss = pd.read_csv(args.in_data_file + '/r.iss.1933pts.txt', sep='\t')[['sample', 'R_ISS', 'Del_17p13.1']]
        df_r2_iss = pd.read_csv(args.in_data_file + '/r2.iss.1933pts.txt', sep='\t')[['sample', 'R2_ISS']]
        df_r_r2_iss = pd.merge(df_r_iss, df_r2_iss, on='sample', how='inner')
        df = pd.merge(df, df_r_r2_iss, on='sample', how='inner')
    else:
        df_r_iss = pd.read_csv(args.in_data_file + '/iss_r-iss_r2-iss_hd6.txt', sep='\t')
        df = pd.merge(df, df_r_iss[['sample', 'R_ISS', 'R2_ISS', 'del17p']], on='sample', how='inner')
        df.rename(columns={'del17p': 'Del_17p13.1'}, inplace=True)

    # GEP70
    if args.get_feats_group == '10':
        df_gep = pd.read_csv(args.in_data_file + '/GEP70.txt', sep='\t')
        df = pd.merge(df, df_gep[['Patient', 'RNA_GEP70_dichotomous']], left_on='sample', right_on='Patient', how='inner')
        df.drop('Patient', axis=1, inplace=True)

    if ext_val is True:
        df['ISS'] = np.select([(df['ISS'] == 'I'), (df['ISS'] == 'II'), (df['ISS'] == 'III')], ['ISS1', 'ISS2', 'ISS3'])
        df["chemo"] = np.select([(df.chemo == 1), (df.chemo == 0)], [0, 0]) # according to the study protocol

    # map clinical vars
    map_gender = {'male': 0, 'female': 1, 'Male':0, 'Female':1}
    x_df = df.replace({'gender': map_gender})
    map_study = {'MMRF': 0, 'MPG': 1, 'Moffit': 2, 'MSKCC_292': 3, 'UAMS': 4, 'HD6': 5}
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

    map_race = {'WHITE': 0, 'BLACK OR AFRICAN AMERICAN': 1, 'OTHER': 2, 'White': 0}
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
    x_df = select_features(x_df, ext_val)
    # print('dataset shape is', x_df.shape)

    return x_df, y_df, df_x_entire


def select_features(df, ext_val=False):
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

    if args.get_feats_group == 'None':
        if (ext_val is not True):
            df_ = df_.drop(['R_ISS', 'R2_ISS'], axis=1)
        elif (ext_val is True):
            df_ = df_.drop(['R_ISS', 'R2_ISS', 'del17p'], axis=1)

    return df_


def raw_data_read_prepare():
    df = pd.read_csv(args.in_data_file + '/PMMM_train_03302023.txt', sep='\t')

    df = df.drop('SCT_line', axis=1)

    df = df.sort_values(by=['sample'])

    x_df, y_df, df_x_entire = data_filter(df, False)

    return x_df, y_df, df_x_entire


def indepn_raw_data_read_prepare():
    '''Raw data for test set'''
    df = pd.read_csv(args.in_data_file + '/PMMM_hd6_03302023.txt', sep='\t')
    df = df.drop('SCT_line', axis=1)

    df = df.sort_values(by=['sample'])

    x_df, y_df, df_x_entire = data_filter(df, True)

    print(x_df.shape, df_x_entire.shape)
    return x_df, y_df, df_x_entire


def cph(train_x, train_y, test_x, val_x_bystate, val_y_bystate):
    train_y = train_y[['binary_event', 'event_time']]
    test_y_bystate = val_y_bystate[['binary_event', 'event_time']]

    cph = CoxnetSurvivalAnalysis(l1_ratio=0.01, fit_baseline_model=True)

    train_x_std_norm, scaler = data_normalize(train_x, 'std_scaler')
    if scaler is not None:
        test_x_std_norm = scaler.transform(test_x)
        test_x_by_state_std_norm = scaler.transform(val_x_bystate)
    else:
        test_x_std_norm = test_x
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


def rsf(train_x, train_y, test_x, val_x_bystate, val_y_bystate, model_param):
    train_y = train_y[['binary_event', 'event_time']]
    test_y_bystate = val_y_bystate[['binary_event', 'event_time']]

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

    surv_val_bystate = model.predict_survival_function(val_x_bystate, return_array=True)
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
        c_ind_val = model.score(val_x_bystate, yval)
    elif args.c_index_type == 'uno':
        try:
            risk_score = model.predict(test_x)
            c_score_uno = concordance_index_ipcw(ytrain, yval, risk_score)
            c_ind_val = c_score_uno[0]
        except Exception:
            c_ind_val = None
            pass

    return surv_val, c_ind_train, c_ind_val, model


def cox_time(train_x, train_y, test_x, val_x_bystate, val_y_bystate, cols_x, model_param, group, df_test):
    np.random.seed(args.seed)
    _ = torch.manual_seed(model_param[2])

    cols_leave = [x for x in cols_x ]
    leave = [(col, None) for col in cols_leave]
    x_mapper = DataFrameMapper(leave)

    train_x = pd.DataFrame(train_x, columns=cols_x)
    train_x = x_mapper.fit_transform(train_x).astype('float32')

    test_x = pd.DataFrame(test_x, columns=cols_x)
    test_x = x_mapper.transform(test_x).astype('float32')

    val_x_bystate = pd.DataFrame(val_x_bystate, columns=cols_x)
    val_x_bystate = x_mapper.transform(val_x_bystate).astype('float32')

    train_x = train_x[:, :len(cols_x)]
    test_x = test_x[:, :len(cols_x)]
    val_x_bystate = val_x_bystate[:, :len(cols_x)]

    train_y = train_y[['binary_event', 'event_time']]
    test_y_bystate = val_y_bystate[['binary_event', 'event_time']]

    labtrans = CoxTime.label_transform()
    get_target = lambda df: (df['event_time'].values, df['binary_event'].values)

    y_train = labtrans.fit_transform(*get_target(train_y))
    y_val = labtrans.transform(*get_target(test_y_bystate))
    val = tt.tuplefy(val_x_bystate, y_val)

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
    surv_test = model.predict_surv_df(test_x)

    df_shap_mean = None
    df_test_shap = None
    c_ind_train = None ; c_ind_val = None

    return surv_test, c_ind_train, c_ind_val, model, df_shap_mean, df_test_shap


def train_survival_(fold, train_x, train_y, val_x, val_y, test_x, df_test, group, model_param):
    # print('train surv', group)

    # aggregated column wise pred probs
    df_pred_probs = pd.DataFrame()

    # feat imps
    df_feat_imp = pd.DataFrame()
    cols_x = list(train_x)

    # impute (loo method)
    #train_x = imputer_knn_loo(train_x)
    #val_x = imputer_knn_loo(val_x)
    #test_x = imputer_knn_loo(test_x)

    # impute
    train_x, imput_model = imputer_knn(train_x, 'train', '')
    val_x, _ = imputer_knn(val_x, 'val', imput_model)
    test_x, _ = imputer_knn(test_x, 'val', imput_model)

    df_train_y = train_y.copy();  df_val_y_bystate = val_y.copy()

    if args.surv_model == 'cph':
        val_surv_df, c_ind_train, c_ind_val, estimator = cph(train_x, df_train_y, test_x, val_x, df_val_y_bystate)

    if args.surv_model == 'rsf':
        val_surv_df, c_ind_train, c_ind_val, estimator = rsf(train_x, df_train_y, test_x, val_x, df_val_y_bystate, model_param)

    if args.surv_model == 'neural_cox_non_prop':
        val_surv_df, c_ind_train, c_ind_val, estimator, mean_shap_values, df_test_shap = \
                                                    cox_time(train_x, df_train_y, test_x, val_x, df_val_y_bystate,
                                                             cols_x, model_param, group, df_test)

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

    df_metrics = pd.DataFrame({'C_index_train': [c_ind_train], 'C_index_val': [c_ind_val], 'Kfolds': [fold],
                               'Fold': [args.kfolds], 'Model': [args.surv_model], 'Group': [group]})

    mean_shap_values, df_test_shap = None, None
    return df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_p1(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list):
    group = 'Progress (P1)'

    tra_x_bystate, tra_y_bystate = crit_prog_p1(train_x, train_y, df_train)
    val_x_bystate, val_y_bystate = crit_prog_p1(val_x, val_y, df_val)

    tra_x_bystate_cp, val_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), val_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        val_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, val_x_bystate_cp, val_y_bystate, test_x_cp, df_test, group, model_param_list[0])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_dec_p1(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list):
    group = 'Progress & deceased (P1)'

    tra_x_bystate, tra_y_bystate = crit_prog_dec_p1(train_x, train_y, df_train)
    val_x_bystate, val_y_bystate = crit_prog_dec_p1(val_x, val_y, df_val)

    tra_x_bystate_cp, val_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), val_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        val_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, val_x_bystate_cp, val_y_bystate, test_x_cp, df_test, group, model_param_list[1])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def move_to_phase2(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list):
    group = 'Move to phase 2'

    tra_x_bystate, tra_y_bystate = crit_move_to_p2(train_x, train_y, df_train)
    val_x_bystate, val_y_bystate = crit_move_to_p2(val_x, val_y, df_val)

    tra_x_bystate_cp, val_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), val_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        val_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, val_x_bystate_cp, val_y_bystate, test_x_cp, df_test, group, model_param_list[2])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_p2(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list):
    group = 'Progress (P2)'

    tra_x_bystate, tra_y_bystate = crit_prog_p2(train_x, train_y, df_train)
    val_x_bystate, val_y_bystate = crit_prog_p2(val_x, val_y, df_val)

    tra_x_bystate_cp, val_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), val_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        val_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, val_x_bystate_cp, val_y_bystate, test_x_cp, df_test, group, model_param_list[3])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def progress_dec_p2(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list):
    group = 'Progress & deceased (P2)'

    tra_x_bystate, tra_y_bystate = crit_prog_dec_p2(train_x, train_y, df_train)
    val_x_bystate, val_y_bystate = crit_prog_dec_p2(val_x, val_y, df_val)

    tra_x_bystate_cp, val_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), val_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        val_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, val_x_bystate_cp, val_y_bystate, test_x_cp, df_test, group, model_param_list[4])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def non_progress_dec_p2(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list):
    group = 'Non-progress & deceased (P2)'

    tra_x_bystate, tra_y_bystate = crit_prog_dec_p2(train_x, train_y, df_train)
    val_x_bystate, val_y_bystate = crit_prog_dec_p2(val_x, val_y, df_val)

    tra_x_bystate_cp, val_x_bystate_cp, test_x_cp = tra_x_bystate.copy(), val_x_bystate.copy(), test_x.copy()

    try:
        tra_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        val_x_bystate_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
        test_x_cp.drop(['SCT_first_line', 'continuos_treat'], axis=1, inplace=True)
    except Exception:
        pass

    df_pred_probs, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = train_survival_(fold, tra_x_bystate_cp,
        tra_y_bystate, val_x_bystate_cp, val_y_bystate, test_x_cp, df_test, group, model_param_list[5])

    return df_pred_probs, group, df_metrics, df_feat_imp, mean_shap_values, df_test_shap


def phase_1(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap):
    # print('phase 11', train_y['os_time'].min(), train_y['os_time'].max())

    # G1 Phase 1 induction (progressors vs non-progressors)
    df_pred_probs_prog_p1, gp_prog_p1, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_p1(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list)
    dict_pred_probs[gp_prog_p1] = df_pred_probs_prog_p1
    dict_feats[gp_prog_p1] = df_feat_imp
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_mean_shap_values[gp_prog_p1] = mean_shap_values
    dict_concat_test_shap[gp_prog_p1] = df_test_shap
    # print('phase 12', train_y['os_time'].min(), train_y['os_time'].max())

    # G2 Phase 1 induction progressors alive vs deceased
    df_pred_probs_prog_dec_p1, gp_prog_dec_p1, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_dec_p1(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list)
    dict_pred_probs[gp_prog_dec_p1] = df_pred_probs_prog_dec_p1
    dict_feats[gp_prog_dec_p1] = df_feat_imp
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_mean_shap_values[gp_prog_dec_p1] = mean_shap_values
    dict_concat_test_shap[gp_prog_dec_p1] = df_test_shap
    # print('phase 13', train_y['os_time'].min(), train_y['os_time'].max())

    # G2 phase 2 vs 1
    df_pred_probs_to_p2, gp_to_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = move_to_phase2(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list)
    dict_pred_probs[gp_to_p2] = df_pred_probs_to_p2
    dict_feats[gp_to_p2] = df_feat_imp
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_mean_shap_values[gp_to_p2] = mean_shap_values
    dict_concat_test_shap[gp_to_p2] = df_test_shap

    return dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap


def phase_2(fold, train_x, train_y, val_x, val_y, test_x, test_y, df_train, df_val, df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap):
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
    df_pred_probs_prog_p2, gp_prog_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_p2(fold, X_p2, Y_p2, val_x, val_y, test_x, test_y, df_p2, df_val, df_test, model_param_list)
    dict_pred_probs[gp_prog_p2] = df_pred_probs_prog_p2
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_feats[gp_prog_p2] = df_feat_imp
    dict_mean_shap_values[gp_prog_p2] = mean_shap_values
    dict_concat_test_shap[gp_prog_p2] = df_test_shap
    #print('phase 22', Y_p2['os_time'].min(), Y_p2['os_time'].max())

    # G5 Phase 2 progress alive vs progress dec
    df_pred_probs_prog_dec_p2, gp_prog_dec_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = progress_dec_p2(fold, X_p2, Y_p2, val_x, val_y, test_x, test_y, df_p2, df_val, df_test, model_param_list)
    dict_pred_probs[gp_prog_dec_p2] = df_pred_probs_prog_dec_p2
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_feats[gp_prog_dec_p2] = df_feat_imp
    dict_mean_shap_values[gp_prog_dec_p2] = mean_shap_values
    dict_concat_test_shap[gp_prog_dec_p2] = df_test_shap
    # print('phase 23', Y_p2['os_time'].min(), Y_p2['os_time'].max())

    # G4 Phase 2 non-prog dec vs rest in p2
    df_pred_probs_non_prog_dec_p2, gp_non_prog_dec_p2, df_metrics, df_feat_imp, mean_shap_values, df_test_shap = non_progress_dec_p2(fold, X_p2, Y_p2, val_x, val_y, test_x, test_y, df_p2, df_val, df_test, model_param_list)
    dict_pred_probs[gp_non_prog_dec_p2] = df_pred_probs_non_prog_dec_p2
    df_concat_metrics = df_concat_metrics.append(df_metrics, ignore_index=True)
    dict_feats[gp_non_prog_dec_p2] = df_feat_imp
    dict_mean_shap_values[gp_non_prog_dec_p2] = mean_shap_values
    dict_concat_test_shap[gp_non_prog_dec_p2] = df_test_shap

    return dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap


def patient_cross_val(model_param_list):
    X, Y, x_df = raw_data_read_prepare()
    X_test, Y_test, x_df_te = indepn_raw_data_read_prepare()

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
    for train_index, val_index in rskf.split(X, x_df_copy['unique_patient_groups']):
        # print(fold, train_index.shape, test_index.shape)

        # ext_val
        train_x, val_x, test_x = X.iloc[train_index].reset_index(drop=True), X.iloc[val_index].reset_index(drop=True), X_test.reset_index(drop=True)
        train_y, val_y, test_y = Y.iloc[train_index].reset_index(drop=True), Y.iloc[val_index].reset_index(drop=True), Y_test.reset_index(drop=True)
        x_df_train = x_df.iloc[train_index].reset_index(drop=True)
        x_df_val = x_df.iloc[val_index].reset_index(drop=True)
        x_df_test = x_df_te.reset_index(drop=True)

        ########## Phase 1 groups ################
        dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap = phase_1(fold,
                                                train_x, train_y, val_x, val_y, test_x, test_y,
                                        x_df_train, x_df_val, x_df_test, model_param_list,
                            dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap)

        ########## Phase 2 groups ################
        dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap = phase_2(fold,
                                    train_x, train_y, val_x, val_y, test_x, test_y, x_df_train, x_df_val,
                            x_df_test, model_param_list, dict_pred_probs, df_concat_metrics, dict_feats, dict_mean_shap_values, dict_concat_test_shap)

        # model probs
        for ind, (k, v) in enumerate(dict_pred_probs.items()):
            if args.get_feats_group == 'None':
                group_id = '100'
            else:
                group_id = args.get_feats_group

            v.to_csv(args.path + '/kfold/model_probs/test_group_' + group_id + '_fold_' + str(fold) + '_' +
                args.surv_model + '_' + k + '.csv', index=False)

        fold += 1

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