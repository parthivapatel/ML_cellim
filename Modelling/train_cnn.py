from numpy.random import seed

seed(42)

import os
import glob
import argparse

import numpy as np

from skimage import io

import sklearn.metrics as sklm
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D, Conv2D
from keras.callbacks import ModelCheckpoint

import warnings

warnings.filterwarnings("ignore")

EPOCH = 200

# c1 max 65535
# c2 max 4671
# c3 max 2442
C_type = 'C3'
max_value = 2442

RADIUS = 32
image_prefix = './t20/2018_02_24_180223_ML_10ngdur_XY'


def load_csv(prefix_folder):

    result = []
    y = []

    for file in glob.glob(prefix_folder):

        with open(file, "r") as f:

            for line in f.readlines():

                line = line.strip().split(",")

                if line[-6] == "20":

                    line.insert(0, file)

                    result.append(line)
                    y.append(int(line[-1]))

    print("There are %s cells in total. " % str(len(result)))
    print()
    return result, y


def load_data(info_time20):

    X_data = []
    y_data = []

    for test_inst in info_time20:

        x_c = int(float(test_inst[1]))
        y_c = int(float(test_inst[2]))

        cell = test_inst[-7]

        file_name = image_prefix + str(cell) + '_' + C_type + '_t20.tif'

        C = io.imread(file_name)

        x_start = x_c - RADIUS
        y_start = y_c - RADIUS

        x_end = x_c + RADIUS
        y_end = y_c + RADIUS

        matrix = C[y_start:y_end, x_start:x_end]

        X_data.append(matrix)
        y_data.append(int(test_inst[-1]))

        if x_start < 0 or y_start < 0:
            raise "Error"

        if x_end > 2048 or y_end > 2044:
            raise "Error"

    X_data = np.array(X_data)
    X_data = np.expand_dims(X_data, axis=3)
    X_data = np.float32(X_data)
    X_data = X_data / max_value

    y_data = np.array(y_data)

    return X_data, y_data


def build_model():

    model = Sequential()

    model.add(Conv2D(8, (5, 5), activation='tanh', input_shape=(64, 64, 1)))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Conv2D(8, (5, 5), activation='tanh'))
    model.add(MaxPooling2D(pool_size=(2, 2)))

    model.add(Flatten())
    model.add(Dense(32, activation='tanh'))
    model.add(Dropout(0.5))

    model.add(Dense(16, activation='tanh'))
    model.add(Dropout(0.5))

    model.add(Dense(1, activation='sigmoid'))

    model.compile(loss='binary_crossentropy',
                  optimizer='adam',
                  metrics=['accuracy'])

    return model


def main(train_index, test_index, csv_data, KFold):

    csv_data = np.array(csv_data)

    directory = "./current_saved_weights/Kfold" + str(KFold)
    file_path = directory + "/weights.hdf5"

    if not os.path.exists(directory):
        os.makedirs(directory)

    X_tr_load = csv_data[train_index]
    X_te_load = csv_data[test_index]

    X_train, y_train = load_data(X_tr_load)
    X_test, y_test = load_data(X_te_load)

    print('Loading complete')
    print()

    checkpoint = ModelCheckpoint(file_path,
                                 monitor='val_loss',
                                 verbose=0,
                                 save_best_only=True)
    callbacks_list = [checkpoint]

    model = build_model()

    history = model.fit(X_train,
                        y_train,
                        epochs=EPOCH,
                        batch_size=32,
                        validation_split=0.3,
                        callbacks=callbacks_list,
                        verbose=0)

    print('Training complete')
    print()

    model.load_weights(file_path)
    scores = np.asarray(model.predict(X_test))

    predictions = []

    for info, score, label in zip(csv_data, scores, y_test):

        key = info[0].split('\\')[-1] + '_' + str(info[-5])
        l = []
        l.append(key)
        l.append(score[0])
        l.append(np.round(score)[0])
        l.append(label)
        predictions.append(l)

    file_name = directory + '/res.csv'
    np.savetxt(file_name, predictions, delimiter=",", fmt="%s")

    print('The accuracy of test is:')
    print(sklm.accuracy_score(y_test, np.round(scores)))
    print()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Parameters of the script')

    parser.add_argument('--k', default=0, help='the kfold')

    parser.add_argument('-csv',
                        default="./180223_ML/*.csv",
                        help='the prefix of csv files')

    parser.add_argument("--random_seed",
                        default=42,
                        help="the randome seed to reproduce the results")

    args = parser.parse_args()

    csv_data, y = load_csv(args.csv)

    skf = StratifiedKFold(n_splits=5,
                          shuffle=True,
                          random_state=args.random_seed)
    KFold = 0

    for train, test in skf.split(csv_data, y):

        main(train, test, csv_data, KFold)

        KFold += 1
