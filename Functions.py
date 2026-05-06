import pandas as pd
import numpy as np


def CreateDifference(dataframe):
    Reference = []
    for i in range(0,16,2):
        Reference.extend([[i for i in range(i,i+7)] for i in range(0,112,7)][i])
    Mutant = []
    for i in range(1,16,2):
        Mutant.extend([[i for i in range(i,i+7)] for i in range(0,112,7)][i])
    number = 0
    something = pd.DataFrame()
    for i in range(0,len(Reference)):
        Difference= dataframe["Feature_"+str(Reference[i])] - dataframe["Feature_"+str(Mutant[i])]
        something.insert(number,"Feature_"+str(number),Difference)
        number += 1
    something.insert(number,"Label",dataframe.iloc[:, 112].values)
    something.insert(number+1,"ID",dataframe.iloc[:, 113].values)
    df = something
    return df


def load_data(files:list,columnsToKeep,Difference = True):
    if type(files) != type([1,2,3]):
        print("files must be a list")
    else:
        dataframes = []
        for file in files:
            df = pd.read_csv(file, sep="\t")
            dataframes.append(df)
        if len(dataframes) > 1:
            df = pd.concat(dataframes, ignore_index=True)
        else:
            df = dataframes[0]
        if Difference:
            df = CreateDifference(df)

        X_train = df.iloc[:, :len(columnsToKeep)].values[:, columnsToKeep]
        y_train = df.iloc[:, len(columnsToKeep)].values

        return [X_train,y_train]


def pad_array(arr, target_length = 30):
    current_length = arr.shape[0]
    if current_length >= target_length:
        return arr  # No padding needed if already at or above target length
    # Calculate padding on both sides
    pad_before = (target_length - current_length) // 2
    pad_after = target_length - current_length - pad_before
    # Apply symmetric padding
    return np.pad(arr, ((pad_before, pad_after), (0, 0)), mode='constant')

def triple_class_predict(modelC, modelD,X_test):
    """
    Combines binary and gain/loss predictions into a single 3-class output.
    Parameters:
        changemodel: Trained Keras model for the prediction of change in binding.
        directionmodel: Trained Keras model for the prediction of type of change in binding.
        X_test: Input data.
    Returns:
        Triple-class predictions: [0 (gained), 1 (lost), 2 (none)]
    """
    binary_preds = modelC.predict(np.abs(X_test))
    binary_class = (binary_preds > 0.5).astype(int).flatten()
    #cols_to_keep = np.concatenate([np.arange(i + 4, i + 7) for i in range(0, 112, 7)])
    #X_filtered = X_test[:, cols_to_keep]
    gain_loss_preds = modelD.predict(X_test)
    gain_loss_class = (gain_loss_preds > 0.5).astype(int).flatten()
    final_preds = np.where(binary_class == 1, 2, gain_loss_class)
    return final_preds


def get_column(infile,column):
    df = pd.read_csv(infile, sep="\t")
    ids = df.iloc[:, column].values
    return ids
