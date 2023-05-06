from sklearn.impute import KNNImputer, SimpleImputer
import sklearn.preprocessing as preprocess
import numpy as np

"""
Helper functions for 
1. Data normalization
2. Imputation - options k-nearest neighbor or mean/median
Please read https://scikit-learn.org/stable/ package documentation for other options
"""


def data_normalize(x_array, normaliz_type):
    """
    :param normaliz_type: data normalize / scalar technique
    :return: preprocessed numpy arrays x,y to feed into ML
    """

    if normaliz_type == 'normalizer':
        scaler = preprocess.Normalizer().fit(x_array)

    # Approach 2 Min-Max
    if normaliz_type == 'min_max':
        # print('Normalization using min-max')
        scaler = preprocess.MinMaxScaler(feature_range=(0, 1))

    if normaliz_type == 'max_abs':
        # print('Max abs scale')
        scaler = preprocess.MaxAbsScaler()

    if normaliz_type == 'std_scaler':
        # print('Standard scaler')
        scaler = preprocess.StandardScaler()

    if normaliz_type == 'power_transformer':
        # print('Power transformer')
        scaler = preprocess.power_transform(method='yeo-johnson')

    if normaliz_type is None:
        # print('No normalization performed')
        return x_array, None

    if x_array.shape[1] != 1:
        x_norm = scaler.fit_transform(x_array)
    else:
        x_norm = scaler.fit_transform(x_array.reshape(-1, 1))
    return x_norm, scaler


def imputer_knn(df, type_, trained_model):
    x_array = df.values
    x_array, _ = data_normalize(x_array, 'min_max')
    if type_ == 'train':
        imputer = KNNImputer(n_neighbors=5, weights='uniform', metric='nan_euclidean')
        imputer.fit(x_array)
        x_imputed = imputer.transform(x_array)
        # print('Missing: %d' % sum(np.isnan(x_imputed).flatten()))
        return x_imputed, imputer
    else:
        x_imputed = trained_model.transform(x_array)
        return x_imputed, _


def imputer_knn_loo(df):
    x_array = df.values
    x_array, _ = data_normalize(x_array, 'min_max')
    imputer = KNNImputer(n_neighbors=5, weights='uniform', metric='nan_euclidean')
    imputer.fit(x_array)
    x_imputed = imputer.transform(x_array)
    # print('Missing: %d' % sum(np.isnan(x_imputed).flatten()))
    return x_imputed


def imputer_simple(df):
    x_array = df.values
    imputer = SimpleImputer(missing_values=np.nan, strategy='median')
    imputer.fit(x_array)
    x_imputed = imputer.transform(x_array)
    # print('Missing: %d' % sum(np.isnan(x_imputed).flatten()))

    return x_imputed

