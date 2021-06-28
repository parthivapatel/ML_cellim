import csv

import warnings
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
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier

warnings.filterwarnings("ignore")

svc_rbf_param_grid = {'C': [0.01, 0.1, 1, 10], 'gamma': [0.001, 0.0001]}

lr_param_grid = {'C': [0.01, 0.1, 1, 10], 'solver': ['liblinear']}

svc_poly_param_grid = {
    'C': [0.01, 0.1, 1, 10],
    'degree': [1, 2, 3],
    'kernel': ['poly']
}


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Parameters of training script")

    parser.add_argument("--cv", default=5, help="the number of folds")

    parser.add_argument("--features",
                        default="v1",
                        help="the features that are used to train the model. v1 is all the features, v2 contains only 3 p65 features. v3 excludes 3 p65 features.")

    parser.add_argument("--random_seed",
                        default=42,
                        help="the randome seed to reproduce the results")

    parser.add_argument("--input_file",
                        default="./ML_2000cells.csv",
                        help="the input file")

    args = parser.parse_args()

    X, y = read_file(args.input_file)

    if args.features == "v1":

        output = "all features"

    elif args.features == "v2":

        X = X[:, 2:5]
        output = "p65 features"

    elif args.features == "v3":

        X = np.delete(X, [2, 3, 4], axis=1)
        output = "none p65 features"

    else:

        raise ("No such options found")

    print(X.shape)

    skf = StratifiedKFold(n_splits=5,
                          shuffle=True,
                          random_state=args.random_seed)

    # Logistic regression
    scores = []

    for train_index, test_index in skf.split(X, y):

        gs = GridSearchCV(LogisticRegression(),
                          param_grid=lr_param_grid,
                          return_train_score=True,
                          scoring='accuracy',
                          cv=skf)

        gs.fit(X[train_index], y[train_index])

        clf = LogisticRegression(**gs.best_params_)

        clf.fit(X[train_index], y[train_index])

        scores.append(accuracy_score(y[test_index],
                                     clf.predict(X[test_index])))

    print("The accuracy of LR using %s: %f" % (output, np.mean(scores)))

    # SVM with rbf kernel
    scores = []

    for train_index, test_index in skf.split(X, y):

        gs = GridSearchCV(SVC(),
                          param_grid=svc_rbf_param_grid,
                          return_train_score=True,
                          scoring='accuracy',
                          cv=skf)

        gs.fit(X[train_index], y[train_index])

        clf = SVC(**gs.best_params_)

        clf.fit(X[train_index], y[train_index])

        scores.append(accuracy_score(y[test_index],
                                     clf.predict(X[test_index])))

    print("The accuracy of SVM-rbf using %s: %f" % (output, np.mean(scores)))

    # SVM with poly kernel
    scores = []

    for train_index, test_index in skf.split(X, y):

        gs = GridSearchCV(SVC(),
                          param_grid=svc_poly_param_grid,
                          return_train_score=True,
                          scoring='accuracy',
                          cv=skf)

        gs.fit(X[train_index], y[train_index])

        clf = SVC(**gs.best_params_)

        clf.fit(X[train_index], y[train_index])

        scores.append(accuracy_score(y[test_index],
                                     clf.predict(X[test_index])))

    print("The accuracy of SVM-poly using %s: %f" % (output, np.mean(scores)))

    # RandomForest
    clf = RandomForestClassifier(n_estimators=100)

    scores = cross_val_score(clf, X, y, cv=skf)

    print("The accuracy of Randomforest using %s: %f" %
          (output, np.mean(scores)))

    # AdaBoost
    clf = AdaBoostClassifier(n_estimators=100)

    scores = cross_val_score(clf, X, y, cv=skf)

    print("The accuracy of AdaBoost using %s: %f" % (output, np.mean(scores)))
