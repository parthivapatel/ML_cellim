import csv

import argparse

import numpy as np
import pandas as pd

from numpy import genfromtxt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV

from sklearn.metrics import accuracy_score

from sklearn.svm import SVC

best_params = {'C': 10, 'gamma': 0.001}


def read_file(file_name):

    data = genfromtxt(file_name, delimiter=',')

    X = data[:, :-1]

    col_mean = np.nanmean(X, axis=0)

    inds = np.where(np.isnan(X))

    X[inds] = np.take(col_mean, inds[1])

    y = data[:, -1]

    scaler = StandardScaler()

    X = scaler.fit_transform(X)

    return X, y


def sample(l, r=1):

    np.random.shuffle(l)

    return l[:int(len(l) * r)]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Parameters of training script")

    parser.add_argument("--cv", default=5, help="the number of folds")

    parser.add_argument("--random_seed",
                        default=42,
                        help="the randome seed to reproduce the results")

    parser.add_argument("--input_file",
                        default="./ML_2000cells.csv",
                        help="the input file")

    parser.add_argument("--num",
                        default=100,
                        help="the number of runs in each fold")

    parser.add_argument("--r",
                        default=0.25,
                        type=float,
                        help="the ratio of training size")

    args = parser.parse_args()

    X, y = read_file(args.input_file)

    skf = StratifiedKFold(n_splits=5,
                          shuffle=True,
                          random_state=args.random_seed)

    accuracy_list = []

    for train_index, test_index in skf.split(X, y):

        fold_list = []

        for _ in range(args.num):

            clf = SVC(C=10, gamma=0.001)

            indices = sample(train_index, args.r)

            clf.fit(X[indices], y[indices])

            y_pred = clf.predict(X[test_index])

            fold_list.append(accuracy_score(y[test_index], y_pred))

        accuracy_list.append(np.mean(fold_list))

    print(np.mean(accuracy_list))