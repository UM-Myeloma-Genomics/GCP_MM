import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Prediction 1.0 Multiple Myeloma workflow')
parser.add_argument('--path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1', help='input dataset path')
parser.add_argument('--legend', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_in/legend_PMMM_2022.xlsx', help='input dataset path')
parser.add_argument('--permut_feat_path', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/permut_feat_imps/aggregate/', help='input dataset path')
parser.add_argument('--permut_feat_path_reranked', type=str, default='/Users/axr2376/Desktop/pred_1_0_paper/data_out/expt_1/permut_feat_imps/aggregate_reranked/', help='input dataset path')
parser.add_argument('--inner_feat_sel', action='store_true', default=False, help='sub features within each group; computationally expensive')
parser.add_argument('--surv_model', type=str, default='neural_cox_non_prop', help='surv model')
parser.add_argument('--iss_feat_groups', action='store_true', default=True, help='iss feature groups')
parser.add_argument('--all_feat_groups', action='store_true', default=True, help='all feature groups')
parser.add_argument('--reverse_order', action='store_true', default=True, help='all feature groups in reverse order')


def create_feat_matrix(df_, feat_gp, top_feat, feat_class='class_generic', reranked=False):
    df_horiz_concat = pd.DataFrame()
    for gp in [('Progress (P1)', top_feat[0]), ('Progress & deceased (P1)', top_feat[1]),
               ('Progress (P2)', top_feat[2]), ('Progress & deceased (P2)', top_feat[3]),
               #('Non-progress & deceased (P2)', top_feat[4]), ('Move to phase 2', top_feat[5])
               ]:
        if reranked:
            df_feat_rank = pd.read_csv(args.permut_feat_path_reranked + gp[0] +
                                        '_' + args.surv_model + '_folds_feat_imp.csv')
        else:
            df_feat_rank = pd.read_csv(args.permut_feat_path + gp[0] +
                                       '_' + args.surv_model + '_folds_feat_imp.csv')

        overall_fea_class = pd.read_excel(args.legend)
        df_feat_rank = pd.merge(df_feat_rank, overall_fea_class[['column_id', 'display_id', 'class_generic', 'class_granular']], on='column_id', how='inner')

        # top n features
        if len(feat_gp.split('|')) == 4:
            split_feat = feat_gp.split('|')
            top_gen_feats = df_feat_rank[(df_feat_rank[feat_class] == split_feat[0])
                                         | (df_feat_rank[feat_class] == split_feat[1])
                                         | (df_feat_rank[feat_class] == split_feat[2])
                                        |  (df_feat_rank[feat_class] == split_feat[3])].iloc[:gp[1]]
        elif len(feat_gp.split('|')) == 1:
            top_gen_feats = df_feat_rank[df_feat_rank[feat_class] == feat_gp].iloc[:gp[1]]

        df_group = pd.DataFrame({gp[0]: top_gen_feats['column_id'].tolist()})

        df_horiz_concat = pd.concat([df_horiz_concat, df_group], axis=1)

    uniqueValues = (df_horiz_concat['Progress (P1)']
                    .append(df_horiz_concat['Progress & deceased (P1)'])
                    .append(df_horiz_concat['Progress (P2)'])
                    .append(df_horiz_concat['Progress & deceased (P2)'])
                    ).unique()

    df_horiz_concat = pd.DataFrame(uniqueValues, columns=['feat_order'])
    df_horiz_concat.dropna(inplace=True)
    df_ = pd.concat([df_, df_horiz_concat], axis=0, ignore_index=True)
    df_ = df_[['feat_order']]

    return df_


def r_iss(df):
    group_name = 'R-ISS'
    ind = 0
    df_feat_combo = pd.DataFrame(columns=['feat_order', 'group_name', 'feat_combo'])
    df_feat_combo.loc[ind, 'feat_combo'] = 'R_ISS'
    df_feat_combo.loc[ind, 'group_name'] = group_name
    df_feat_combo.drop('feat_order', axis=1, inplace=True)
    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def r2_iss(df):
    group_name = 'R2-ISS'
    ind = 0
    df_feat_combo = pd.DataFrame(columns=['feat_order', 'group_name', 'feat_combo'])
    df_feat_combo.loc[ind, 'feat_combo'] = 'R2_ISS'
    df_feat_combo.loc[ind, 'group_name'] = group_name
    df_feat_combo.drop('feat_order', axis=1, inplace=True)
    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss(df):
    group_name = 'ISS'
    ind = 0
    df_feat_combo = pd.DataFrame(columns=['feat_order', 'group_name', 'feat_combo'])
    df_feat_combo.loc[ind, 'feat_combo'] = 'ISS'
    df_feat_combo.loc[ind, 'group_name'] = group_name
    df_feat_combo.drop('feat_order', axis=1, inplace=True)
    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def gep70(df):
    group_name = 'RNA_GEP70_dichotomous'
    ind = 0
    df_feat_combo = pd.DataFrame(columns=['feat_order', 'group_name', 'feat_combo'])
    df_feat_combo.loc[ind, 'feat_combo'] = 'RNA_GEP70_dichotomous'
    df_feat_combo.loc[ind, 'group_name'] = group_name
    df_feat_combo.drop('feat_order', axis=1, inplace=True)
    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_clinical(df):
    group_name = 'ISS-Clinical'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_clinical_transloc(df):
    group_name = 'ISS-Clinical-RecTransloc'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4*[155], 'class_granular')

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_clinical_transloc_sct(df):
    group_name = 'ISS-Clinical-RecTransloc-SCT'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4*[155], 'class_granular')
    df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_clinical_transloc_sct_conTreat(df):
    group_name = 'ISS-Clinical-RecTransloc-SCT-ContTreat'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4*[155], 'class_granular')
    df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_clinical_transloc_sct_conTreat_treatment(df):
    group_name = 'ISS-Clinical-RecTransloc-SCT-ContTreat-Treatment'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4*[155], 'class_granular')
    df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_clinical_transloc_sct_conTreat_treatment_gen(df, perf_type):
    group_name = 'ISS-Clinical-RecTransloc-SCT-ContTreat-Treatment-TopGen'
    df_feat_stack = pd.DataFrame()

    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4*[155])

    df_uni = pd.read_csv(args.path + '/plots/df_feat_uni_p_val_05.csv')

    if perf_type == 'top5':
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations', [5, 2, 5, 5],
                                           'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])
    elif perf_type == 'top10':
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations',
                                           [10, 2, 10, 10], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])
    elif perf_type == 'uni':
        print('just univariate')
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    elif perf_type == 'uni_p_01':
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations',
                                           [5, 2, 5, 5], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    elif perf_type == 'uni_p_01_only':
        print('just univariate')
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    elif perf_type == 'uni_p_01_top5':
        print('just univariate')
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations',
                                           [1, 1, 3, 2], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    elif perf_type == 'uni_p_01_top10':
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations',
                                           [9, 9, 10, 10], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    elif perf_type == 'uni_p_01_top5_from_p_05':
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations|Translocation',
                                           [2, 2, 3, 4], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_transloc(df):
    group_name = 'ISS-RecTransloc'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4*[155], 'class_granular')

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_sct(df):
    group_name = 'ISS-SCT'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_cont_treat(df):
    group_name = 'ISS-Continuous treatment'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_therapy(df):
    group_name = 'ISS-Treatment'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])
    df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4*[155])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def iss_gen(df, perf_type):
    group_name = 'ISS-TopGen'
    df_feat_stack = pd.DataFrame()
    df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4*[155])

    df_uni = pd.read_csv(args.path + '/plots/df_feat_uni_p_val_05.csv')

    if perf_type == 'uni_p_01_top5_from_p_05':
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations|Translocation',
                                           [2, 1, 3, 4], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def reverse_order_feats(df, perf_type, conc=None):
    df_feat_stack = pd.DataFrame()
    group_name = 'TopGen'
    if conc == 'therapy':
        group_name = 'TopGen+Therapy'
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4 * [155])
    elif conc == 'maintainence':
        group_name = 'TopGen+Therapy+Maintainence'
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4 * [155])
    elif conc == 'sct':
        group_name = 'TopGen+Therapy+Maintainence+SCT'
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4 * [155])
    elif conc == 'transloc':
        group_name = 'TopGen+Therapy+Maintainence+SCT+Transloc'
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4 * [155], 'class_granular')
    elif conc == 'clinical':
        group_name = 'TopGen+Therapy+Maintainence+SCT+Transloc+Clinical'
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4 * [155], 'class_granular')
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4 * [155])
    elif conc == 'iss':
        group_name = 'TopGen+Therapy+Maintainence+SCT+Transloc+Clinical+ISS'
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Continuous treatment', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'SCT', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Translocation', 4 * [155], 'class_granular')
        df_feat_stack = create_feat_matrix(df_feat_stack, 'Clinical', 4 * [155])
        df_feat_stack = create_feat_matrix(df_feat_stack, 'ISS', 4 * [155])

    df_uni = pd.read_csv(args.path + '/plots/df_feat_uni_p_val_05.csv')

    if perf_type == 'uni_p_01_only':
        print('just univariate')
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    elif perf_type == 'uni_p_01_top5_from_p_05':
        df_uni = df_uni[(df_uni['PFS.pvalue'] < 0.01) & (df_uni['OS.pvalue'] < 0.01)]
        df_feat_stack = create_feat_matrix(df_feat_stack, 'CNV|Mutations|CNV_Mutations|Translocation',
                                           [2, 2, 3, 4], 'class_granular', True)
        df_feat_stack = pd.DataFrame(df_feat_stack['feat_order'].append(df_uni['Genomic.Feature'], ignore_index=True),
                     columns=['feat_order'])

    feats = df_feat_stack['feat_order'].to_list()
    cat_elem = []
    for ind, item in enumerate(feats):
        if args.inner_feat_sel and ind < len(feats)-1:
                cat_elem.append(feats[:ind+1])
        elif ind == len(feats)-1:
            cat_elem.append(feats[:ind + 1])

    cat_elem = [list(set(cat_elem[0]))]

    cat_elem = [' '.join(item for item in i) for i in cat_elem]

    df_feat_combo = pd.DataFrame(cat_elem, columns=['feat_combo'])
    df_feat_combo['group_name'] = group_name

    df = pd.concat([df, df_feat_combo], ignore_index=True)
    return df


def main():
    global args
    args = parser.parse_args()

    df = pd.DataFrame(columns=['feat_combo', 'group_name'])
    if args.iss_feat_groups:
        df = r_iss(df)
        df = r2_iss(df)
        df = iss(df)

    if args.all_feat_groups:
        df = iss_clinical(df)
        df = iss_clinical_transloc(df)
        df = iss_clinical_transloc_sct(df)
        df = iss_clinical_transloc_sct_conTreat(df)
        df = iss_clinical_transloc_sct_conTreat_treatment(df)
        df = iss_clinical_transloc_sct_conTreat_treatment_gen(df, 'uni_p_01_only')
        df = iss_clinical_transloc_sct_conTreat_treatment_gen(df, 'uni_p_01_top5_from_p_05')
        ###
        df = gep70(df)
        df = iss_sct(df)
        df = iss_cont_treat(df)
        df = iss_transloc(df)
        df = iss_therapy(df)
        df = iss_gen(df, 'uni_p_01_top5_from_p_05')

    if args.reverse_order:
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05')
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05', 'therapy')
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05', 'maintainence')
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05', 'sct')
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05', 'transloc')
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05', 'clinical')
        df = reverse_order_feats(df, 'uni_p_01_top5_from_p_05', 'iss')

    df = df[~( (df['feat_combo'] == 'ISS') & (df['group_name'] != 'ISS') )].reset_index(drop=True)
    df['group_id'] = df.index

    # main file
    df.to_csv(args.path + '/feat_matrix/main_keeper.csv', index=False)

    return None


if __name__ == '__main__':
    main()