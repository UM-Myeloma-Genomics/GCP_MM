import numpy as np

"""
Helper functions for annotations of events and event times
1. Progression (P1)
2. Progress & dec. (P1)
3. Move to phase 2
4. Progress (P2)
5. Progress & dec (P2)
6. Non-progress & dec (P2)
 
Please read paper/supplementary for more details 
"""


def crit_prog_p1(t_x, t_y, df):
    y_p1 = t_y.copy()
    y_p1['event_time'] = np.where(df['pfs_time'] > df['duration'], df['duration'], df['pfs_time'])

    subset_df_class1 = df[(df['phases'] == 'phase_1') & (df['pfs_code'] == 1)]
    y_p1.loc[subset_df_class1.index.tolist(), 'binary_event'] = 1

    subset_df_class0 = df[~df.index.isin(subset_df_class1.index)]
    y_p1.loc[subset_df_class0.index.tolist(), 'binary_event'] = 0
    return t_x, y_p1


def crit_prog_dec_p1(t_x, t_y, df):
    y_prog_dec_1 = t_y.copy()

    df_X_prog_p1 = df[(df['phases'] == 'phase_1') & (df['pfs_code'] == 1)]
    t_x = t_x[t_x.index.isin(df_X_prog_p1.index)]
    y_prog_dec_1 = y_prog_dec_1[y_prog_dec_1.index.isin(df_X_prog_p1.index)]

    y_prog_dec_1['os_time'] = y_prog_dec_1['os_time'] - df['pfs_time']  # reset

    subset_df_class1 = df_X_prog_p1[(df_X_prog_p1['os_code'] == 0)]
    y_prog_dec_1.loc[subset_df_class1.index.tolist(), 'binary_event'] = 0
    subset_df_class0 = df_X_prog_p1[(df_X_prog_p1['os_code'] == 1)]
    y_prog_dec_1.loc[subset_df_class0.index.tolist(), 'binary_event'] = 1
    y_prog_dec_1['event_time'] = y_prog_dec_1['os_time']

    return t_x, y_prog_dec_1


def crit_move_to_p2(t_x, t_y, df):
    y_p1_to_2 = t_y.copy()

    df['duration'] = np.where(df['duration'] >= 365, 365, df['duration'])

    subset_df_class1 = df[(df['phases'] == 'phase_2')]
    y_p1_to_2.loc[subset_df_class1.index.tolist(), 'binary_event'] = 1
    subset_df_class0 = df[(df['phases'] == 'phase_1')]
    y_p1_to_2.loc[subset_df_class0.index.tolist(), 'binary_event'] = 0

    y_p1_to_2['event_time'] = np.where(df['pfs_time'] > df['duration'], df['duration'], df['pfs_time'])

    return t_x, y_p1_to_2


def crit_prog_p2(t_x, t_y, df):
    y_p2 = t_y.copy()
    df_cp = df.copy()

    subset_df_class1 = df_cp[(df_cp['pfs_code'] == 1)]
    y_p2.loc[subset_df_class1.index.tolist(), 'binary_event'] = 1
    subset_df_class0 = df_cp[(df_cp['pfs_code'] == 0)]
    y_p2.loc[subset_df_class0.index.tolist(), 'binary_event'] = 0
    y_p2['event_time'] = y_p2['pfs_time']

    return t_x, y_p2


def crit_prog_dec_p2(t_x, t_y, df):
    y_prog_dec_p2 = t_y.copy()

    df_X_prog_p2 = df[(df['phases'] == 'phase_2') & (df['pfs_code'] == 1)]
    t_x = t_x[t_x.index.isin(df_X_prog_p2.index)]
    y_prog_dec_p2 = y_prog_dec_p2[y_prog_dec_p2.index.isin(df_X_prog_p2.index)]

    subset_df_class1 = df_X_prog_p2[(df_X_prog_p2['pfs_code'] == 1) & (df_X_prog_p2['os_code'] == 0)]
    y_prog_dec_p2.loc[subset_df_class1.index.tolist(), 'binary_event'] = 0
    subset_df_class0 = df_X_prog_p2[(df_X_prog_p2['pfs_code'] == 1) & (df_X_prog_p2['os_code'] == 1)]
    y_prog_dec_p2.loc[subset_df_class0.index.tolist(), 'binary_event'] = 1

    # reset
    y_prog_dec_p2['event_time'] = y_prog_dec_p2['os_time'] - df_X_prog_p2['pfs_time'] + df_X_prog_p2['duration']

    return t_x, y_prog_dec_p2


def crit_non_prog_dec_p2(t_x, t_y, df):
    subset_df_class1 = df[(df['pfs_code'] == 0) & (df['os_code'] == 1)]
    t_y.loc[subset_df_class1.index.tolist(), 'binary_event'] = 1
    subset_df_class0 = df[~df.index.isin(subset_df_class1.index)]
    t_y.loc[subset_df_class0.index.tolist(), 'binary_event'] = 0
    censor_prog = df[(df['pfs_code'] == 1)]
    t_y.loc[censor_prog.index.tolist(), 'os_time'] = censor_prog['pfs_time']
    t_y['event_time'] = t_y['os_time']
    return t_x, t_y

