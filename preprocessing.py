from sklearn.preprocessing import OneHotEncoder, FunctionTransformer, MinMaxScaler
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.compose import make_column_transformer
from sklearn.pipeline import Pipeline
from sklearn import set_config
import pandas as pd
import numpy as np

set_config(transform_output="pandas")


def make_onehot(onehot_columns: list):
    onehot = make_column_transformer(
        (
            OneHotEncoder(handle_unknown="ignore", sparse_output=False, dtype=np.int8),
            onehot_columns,
        ),
        remainder="passthrough",
        verbose_feature_names_out=False,
    )
    return onehot


def make_imputer(impute_columns: list):
    imputer = make_column_transformer(
        (
            KNNImputer(weights="distance", n_neighbors=5, add_indicator=True),
            impute_columns,
        ),
        remainder="passthrough",
        verbose_feature_names_out=False,
    )
    return imputer


def make_scaler(scaler_columns: list):
    scaler = make_column_transformer(
        (
            MinMaxScaler(),
            scaler_columns,
        ),
        remainder="passthrough",
        verbose_feature_names_out=False,
    )
    return scaler
    

def Create_Pipeline(cols_to_encode, cols_to_scale, imputation=False, cols_to_impute=None):
    onehot = make_onehot(cols_to_encode)
    scale = make_scaler(cols_to_scale)
    if imputation:
        impute = make_imputer(cols_to_impute)
        preprocess_pipeline = Pipeline(
            [
            ("onehot", onehot),
            ("impute", impute),
            ("scale", scale)
            ]
        )
    else:
        preprocess_pipeline = Pipeline(
            [
            ("onehot", onehot),
            ("scale", scale)
            ]
        )

    return preprocess_pipeline


if __name__ == "__main__":
    cols_to_encode = ['REF', 'ALT', 'GENE_SYMBOL', 'TYPE', 'VARIANT_TYPE', 
                      'MOST_SEVERE_CONSEQUENCE', 'IMPACT', 'EXON_INTRON_TYPE', 'CLINPRED_PRED']
    cols_to_impute = ['CLINPRED_RANKSCORE', 'POLYPHEN2_HDIV_RANKSCORE',  'SIFT_CONVERTED_RANKSCORE',
                      'SIFT4G_CONVERTED_RANKSCORE',  'MUTATIONASSESSOR_RANKSCORE', 'MUTATIONTASTER_CONVERTED_RANKSCORE']
    cols_to_scale = ['STRAND', 'VARIANT_OCCURRENCES', 'EXON_INTRON_N', 'DOMAINS_COUNT', 'PUBMED_COUNT']

    pipeline = Create_Pipeline(cols_to_encode, cols_to_scale, imputation=True, cols_to_impute=cols_to_impute)
    df = pd.read_csv("data/ML_BRCA.csv")
    df = pipeline.fit_transform(df)
    print(df.head())
    print(df.shape)
    print(df.dtypes)
    print(df.isna().sum())

    df.to_csv("data/ML_BRCA_preprocessed.csv", index=False)