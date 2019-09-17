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
import numpy as np
from sklearn.ensemble import RandomForestClassifier


class ForestTrainer:

    def __init__(self, training_file_path, k = 50, stop_parameter = 0.01):

        self.__k = k
        self.__stop_parameter = stop_parameter

        data_holder = self.__load_observation_data(training_file_path)

        self.__test_data = data_holder['testing_data']
        self.__train_data = data_holder['training_data']
        self.__validation_data = data_holder['validation_data']

        self.__clf = self.__train_forest()

    def get_forest(self):

        return self.__clf

    def __train_forest(self):

        ## Use the number of true calls to set out our even class balance:

        ## Now run the random forest until we hit our stop parameter
        first_forest = self.__run_random_forest()

        perf_diff = 1.00
        last_perf = first_forest['accuracy']
        current_forest = first_forest['forest']

        performance_tracker = dict()
        current_iteration = 1
        performance_tracker[current_iteration] = last_perf

        while perf_diff > self.__stop_parameter and self.__validation_data.shape[0] > self.__k:

            current_iteration = current_iteration + 1
            recursive_forest = self.__run_random_forest()

            perf_diff = abs(last_perf - recursive_forest['accuracy'])
            last_perf = recursive_forest['accuracy']

            performance_tracker[current_iteration] = last_perf

            current_forest = recursive_forest['forest']

        return current_forest

    def __run_random_forest(self):

        value_columns = self.__train_data.columns[2:19]
        train_values = self.__train_data[value_columns]
        train_outcome = self.__train_data["annotation"]

        clf = RandomForestClassifier(n_estimators=500)
        clf.fit(train_values, train_outcome)

        validation_values = self.__validation_data[value_columns]

        validation_predictions = clf.predict_proba(validation_values)

        new_validation_data = self.__validation_data.copy()

        new_validation_data["prob_N"] = validation_predictions.T[0]
        new_validation_data["prob_Y"] = validation_predictions.T[1]

        new_validation_data["diff_probs"] = abs(new_validation_data["prob_N"] - new_validation_data["prob_Y"])

        add_to_train = new_validation_data.sort_values(by='diff_probs')[0:self.__k].copy()
        add_to_train = add_to_train.drop(["prob_N", "prob_Y", "diff_probs"], axis=1)

        self.__train_data = pd.concat([self.__train_data, add_to_train])

        self.__validation_data = new_validation_data[~new_validation_data.index.isin(add_to_train.index)]

        test_values = self.__test_data[value_columns]
        test_predictions = clf.predict_proba(test_values)
        self.__test_data["prob_Y"] = test_predictions.T[1]

        ## Calculate Current Accuracy
        true_positives = self.__test_data[self.__test_data.annotation == True]
        num_tp = float(true_positives[true_positives.prob_Y >= 0.5].shape[0])
        true_negatives = self.__test_data[self.__test_data.annotation == False]
        num_tn = float(true_negatives[true_negatives.prob_Y < 0.5].shape[0])

        num_p = float(self.__test_data[self.__test_data.annotation == True].shape[0])
        num_n = float(self.__test_data[self.__test_data.annotation == False].shape[0])

        accuracy = (num_tp + num_tn) / (num_p + num_n)

        returnable = dict()

        returnable['accuracy'] = accuracy
        returnable['forest'] = clf

        return returnable

    def __load_observation_data(self, training_file_path):

        observation_data = pd.read_csv(training_file_path, sep="\t")

        ## Check possible values of annotation and modify if necessary:
        poss_values = observation_data.annotation.unique()
        if observation_data.dtypes['annotation'] is np.dtype(np.bool):
            print ("Evaluating training data as raw (True/False) boolean annotations...\n")

        elif observation_data.dtypes['annotation'] is np.dtype(np.int64):
            ## Check if just 0 and 1 and convert to boolean
            if 0 in poss_values and 1 in poss_values and poss_values.size == 2:
                print ("Evaluating training data as converted binary (0/1) boolean annotations...\n")
                observation_data['annotation'] = observation_data['annotation'].astype('bool')
            else:
                print ("Annotation column either has only 0 or 1, or has numbers other than 0... Exiting...\n")
                exit(1)

        else:
            ## check for only T/F in annotation
            if "T" in poss_values and "F" in poss_values and poss_values.size == 2:
                print ("Evaluating training data as shortened (T/F) boolean annotations...\n")
                observation_data['annotation'] = observation_data["annotation"] == "T"

            else:
                print ("Annotation column did not contain possible expected boolean values... exiting...\n")
                exit(1)

        num_true = observation_data[observation_data.annotation == True].shape[0]

        observation_data_even = pd.concat([observation_data[observation_data.annotation == True],
                                           observation_data[observation_data.annotation == False].sample(
                                               n=num_true, replace=False)])

        train_data = observation_data_even.sample(n=self.__k, replace=False)
        used_rows = list(train_data.index)

        test_data = observation_data_even[~observation_data_even.index.isin(used_rows)].sample(n=self.__k, replace=False)
        used_rows.extend(list(test_data.index))

        validation_data = observation_data_even[~observation_data_even.index.isin(used_rows)]

        data_holder = dict()
        data_holder['training_data'] = train_data
        data_holder['testing_data'] = test_data
        data_holder['validation_data'] = validation_data

        return data_holder