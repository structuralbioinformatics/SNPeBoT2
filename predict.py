import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, Model, regularizers
from tensorflow.keras.optimizers import AdamW
import keras_tuner as kt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
import random, os
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from tensorflow.keras.layers import Attention, Reshape
import tensorflow.keras.backend as K
from sklearn.metrics import classification_report
from sklearn.utils.class_weight import compute_class_weight
import pickle
from keras.models import Sequential
from keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
from collections import Counter
import Functions as Fc
import sys

# Set random seeds for reproducibility
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)
os.environ['PYTHONHASHSEED'] = str(SEED)

Input=sys.argv[1]
Output=sys.argv[2]

# Extract Features and Labels
cols_to_keep = np.array([i for i in range(0, 56)]) # KEEP ALL DATA

X_test,y_test = Fc.load_data([Input],cols_to_keep)#Resources/Balanced_Test.tsv

Ids = Fc.get_column(Input,113)


CNmodel = keras.models.load_model(f"Models/BenchPlus_ChangevsNone_General.keras")
GLmodel = keras.models.load_model(f"Models/BenchPlus_GainvsLoss_General.keras")


output = Fc.triple_class_predict(CNmodel,GLmodel,X_test)
#print(Ids)
#print(output)

translateOutput = {0:"gain",1:"loss",2:"no-change"}

wt=open(Output,"w")

for i in range(0,len(Ids)):
    addition = Ids[i] + "\t" + translateOutput[output[i]] + "\n"
    wt.write(addition)

wt.close()




