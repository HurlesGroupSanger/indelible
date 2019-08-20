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


def score_file(forest_path, testing_file_path):

    df = pd.read_csv(testing_file_path, sep="\t")
    values = df.loc[:, df.columns.difference(["chrom", "position", "seq_longest", "pct_double_split"])]
    clf = loadForest(forest_path)


    predicted_class = clf.predict(values)
    df["predicted"] = predicted_class
    df["prob_N"] = clf.predict_proba(values).T[0]
    df["prob_Y"] = clf.predict_proba(values).T[1]
    return df


def saveForest(clf, output_path):
    with bz2.BZ2File(output_path, 'wb') as fid:
        pickle.dump(clf, fid)
    return True


def loadForest(path):
    with bz2.BZ2File(path, 'rb') as fid:
        clf = pickle.load(fid)
    return clf


def score_positions(input_path, output_path, config):
    df = score_file(config['random_forest_model'], input_path)
    df.to_csv(output_path, sep="\t", index=False)


def train(input_path, output_path):
    clf = trainForest(input_path)
    saveForest(clf, output_path)
