"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute

    Loads the Random Forest model and scores all the variants in the file

    Parameters
    ----------

    1) File to be scored
	2) Model file

    Returns
    -------

    1) File with appended scored

"""
import sys
import timeit

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import pickle
import bz2


def trainForest(training_file_path):
    df = pd.read_csv(training_file_path, sep="\t")
    values = df.ix[:, df.columns - ["chrom", "position", "annotation", "sample", "seq_longest", "pct_double_split"]]
    outcome = df["annotation"]
    clf = RandomForestClassifier(n_estimators=10000)
    clf.fit(values, outcome)
    return clf


def score_dataframe(clf, df):

    values = df.loc[:, df.columns.difference(["chrom", "position", "seq_longest", "pct_double_split"])]

    df["prob_N"] = clf.predict_proba(values).T[0]
    df["prob_Y"] = clf.predict_proba(values).T[1]
    return df


def calculate_prediction_column(row):

    if row.prob_N > row.prob_Y:
        return "N"
    else:
        return "Y"

def saveForest(clf, output_path):
    with bz2.BZ2File(output_path, 'wb') as fid:
        pickle.dump(clf, fid)
    return True


def loadForest(path):
    with bz2.BZ2File(path, 'rb') as fid:
        clf = pickle.load(fid)
    return clf


def score_positions(input_path, output_path, config):

    clf = loadForest(config['random_forest_model'])

    df = pd.read_csv(input_path, sep="\t")
    df_length = len(df.index)
    df_final = pd.DataFrame()

    if df_length <= 20000:
        df_final = score_dataframe(clf, df)
    else:
        # Need to run df as chunks through the file:
        for i in range(0,df_length,20000):
            if i + 20000 > df_length:
                df_chunk = df[i:df_length]
                df_chunk = df_chunk.copy()
            else:
                df_chunk = df[i:i+20000]
                df_chunk = df_chunk.copy()

            df_final = df_final.append(score_dataframe(clf,df_chunk), ignore_index=True)

    df_final["predicted"] = df_final.apply(lambda row: calculate_prediction_column(row), axis=1)
    df_final.to_csv(output_path, sep="\t", index=False)


def train(input_path, output_path):
    clf = trainForest(input_path)
    saveForest(clf, output_path)
