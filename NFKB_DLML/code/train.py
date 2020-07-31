from numpy.random import seed
from tensorflow import set_random_seed

seed(42)
set_random_seed(42)

import os 
import glob
import pickle 
from collections import Counter

import pandas as pd
import numpy as np

import skimage
from skimage import io

import sklearn.metrics as sklm
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.linear_model import LogisticRegression

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D, Conv2D
from keras.utils import np_utils
from keras.callbacks import ModelCheckpoint

import warnings
warnings.filterwarnings("ignore")

# KFold number should be from 0 to 9
KFold = 0 

EPOCH = 200
max_value = 2442
C_type = 'C3'
RADIUS = 32
image_prefix = './image_t20/2018_02_24_180223_ML_10ngdur_XY'
file_path="./current_saved_weights/Kfold" + str(KFold) + "/weights-improvement-{epoch:02d}-{val_acc:.2f}.hdf5"


def load_data(info_time20):
    
    X_data = []
    y_data = []

    for test_inst in info_time20:

        x_c = int(test_inst[1])
        y_c = int(test_inst[2])

        cell = test_inst[-7]

        file_name = image_prefix + str(cell) + '_' + C_type + '_t20.tif'

        C = io.imread(file_name)

        x_start = x_c - RADIUS
        y_start = y_c - RADIUS

        x_end = x_c + RADIUS
        y_end = y_c + RADIUS 

        matrix = C[y_start:y_end, x_start:x_end]

        X_data.append(matrix)
        y_data.append(test_inst[-1])

        if x_start < 0 or y_start < 0:
            raise "Error"

        if x_end > 2048 or y_end > 2044:
            raise "Error"
            
    X_data = np.array(X_data)       
    X_data = np.expand_dims(X_data, axis=3)
    X_data = np.float32(X_data)
    X_data = X_data/max_value
    
    y_data = np.array(y_data)
    
    return X_data, y_data


def build_model():
    
    model = Sequential()
    
    model.add(Conv2D(8, (5, 5), activation='tanh', input_shape=(64,64,1)))
    model.add(MaxPooling2D(pool_size=(2,2)))
    model.add(Conv2D(8, (5, 5), activation='tanh'))
    model.add(MaxPooling2D(pool_size=(2,2)))

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


def main():

    directory = "./current_saved_weights/Kfold" + str(KFold) 

    if not os.path.exists(directory):
        os.makedirs(directory)

    with open('./pickle/index.pickle', 'rb') as f:
        indices = pickle.load(f)
        
    train, valid, test = indices[KFold]

    with open('./pickle/data.pickle', 'rb') as f:
        data = pickle.load(f)
    
    X = data['X']
    
    X_tr_load = X[train]
    X_va_load = X[valid]
    X_te_load = X[test]

    info_tr_time20 = X_tr_load[:,4,:]
    info_va_time20 = X_va_load[:,4,:]
    info_te_time20 = X_te_load[:,4,:]

    X_train, y_train = load_data(info_tr_time20)
    X_valid, y_valid = load_data(info_va_time20)
    X_test,  y_test  = load_data(info_te_time20)

    
    ratio_positive = sum(y_train)/len(y_train)
    ratio_negative = 1 - ratio_positive

    weighted_dict = {}
    weighted_dict ={1:ratio_negative, 0:ratio_positive}

    checkpoint = ModelCheckpoint(file_path, monitor='val_acc', verbose=0, save_best_only=False)
    callbacks_list = [checkpoint]

    model = build_model()

    history = model.fit(X_train, y_train, epochs=EPOCH, 
                        batch_size=32, 
                        class_weight = weighted_dict,
                        validation_data=(X_valid, y_valid), 
                        callbacks=callbacks_list,
                        verbose=1)

    idx = np.argmax(history.history['val_acc'])

    print("The best epoch number is %d (starting from 1)" %(idx +1))


if __name__ == "__main__":

    main()

