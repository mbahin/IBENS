import argparse
import matplotlib.pyplot as plt
#matplotlib.use("Agg")
import sys
import numpy as np
import keras
from keras.layers import Dense, Flatten, InputLayer, Conv1D, MaxPooling1D
from keras.models import Sequential
import pandas as pd

# Genes: ["BMP4", "Bactine2", "Crabp2a", "DeltaC", "DeltaD", "FGF8a", "Her7", "Hoxa1a", "Hoxa5a", "Hoxb1b", "Hoxd4a", "Notch1a", "RPL13a1", "RPL13a2", "Rarab"]
# On research plate:
## PS: regular doses
## PC: variable doses (less reliable)
# Not reliable: Notch1a/PC, BMP4/PC & Crabp2a
# PS: Notch1a

def convert(y):
    my = np.mean(y)
    sy = np.std(y)
    return (y - my) / sy, my, sy


def back_convert(y, my, sy):
    return y * sy + my


def read_file(df, test_dataset=None, validation=False, exclude=None):
    ''' Function to load the data from the input file and dividing the dataset into train/test if necessary.'''
    # Loading the data
    gene = np.array(df["Gene"])
    condition = np.array(df["Condition"]).astype(str)
    wellID = np.array(df["WellID"])
    y = np.array(df["Theoretical_N0_log"])
    x = np.array(df.loc[:, df.columns.str.startswith("F")])

    # Excluding gene/conditions if required
    if exclude != [()]:
        idx_del = np.logical_or.reduce([((gene == g) & np.char.startswith(condition, c)) for (g, c) in exclude])
        gene = gene[~idx_del]
        condition = condition[~idx_del]
        wellID = wellID[~idx_del]
        y = y[~idx_del]
        x = x[~idx_del]

    # Changing the range of values to be low and around 0
    y, my, sy = convert(y)

    if validation:
        # If this is the validation step, there is no train/test, the full dataset is returned
        return wellID, (gene, condition, np.expand_dims(x, axis=2), y), my, sy
    else:
        # Getting the test indexes
        idx_test = np.logical_or.reduce([((gene == g) & np.char.startswith(condition, c)) for (g, c) in test_dataset])
        # Creating the train dataset
        gene_train = gene[~idx_test]
        condition_train = condition[~idx_test]
        x_train = np.expand_dims(x[~idx_test], axis=2)
        y_train = y[~idx_test]
        # Creating the test dataset
        gene_test = gene[idx_test]
        condition_test = condition[idx_test]
        x_test = np.expand_dims(x[idx_test], axis=2)
        y_test = y[idx_test]
        # Returning train and test datasets
        return wellID, (gene_train, condition_train, x_train, y_train), (gene_test, condition_test, x_test, y_test), my, sy


#def display_plot(model, x, y, my, sy, gene, condition):
def display_plot(model, x, y, my, sy, gene):
    ''' Function to display the theoretical vs the predicted. '''
    # Setting the colors
    hsv = plt.get_cmap("hsv")
    color = hsv(np.linspace(0, 1.0, np.unique(gene).shape[0] + 1))
    # Converting back the data to the original range
    real = np.log(back_convert(y, my, sy))
    pred = np.log(back_convert(model.predict(x).T[0], my, sy))
    # Plotting each gene results
    for i, g in enumerate(np.unique(gene)):
        idx = np.where(gene == g)
        #plt.scatter(real[idx], pred[idx], label=g, c=color[i], edgecolors=["k" if b else color[i] for b in np.char.startswith(condition[idx], "PC")])  # Circling the test PC homologue
        plt.scatter(real[idx], pred[idx], label=g, c=color[i])
    plt.plot([-3, 3], [-3, 3], "-")
    plt.legend(loc=0)
    plt.grid("on")
    plt.axis("equal")
    plt.show()


# Main
if __name__ == "__main__":

    # Setting the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config1", action="store_true", dest="config1", help="Plate2 is used for training/testing and plate1 for validation.")
    parser.add_argument("--config2", action="store_true", dest="config2", help="Plate1 is used for training/testing and plate2 for validation.")
    options = parser.parse_args()

    # Setting variables
    cycles_col = ["F" + str(i).zfill(2) for i in range(40)]

    # Setting gene/condition lists
    test_gene_condition = [("FGF8a", "PS"), ("Rarab", "PS"), ("Hoxd4a", "PS")]
    #test_gene_condition = [("FGF8a", "PS"), ("Her7", "PC"), ("Crabp2a", "PS")]  # Only found test set where val_loss > loss from the beginning for config1
    # Loading plate1 dataset
    QV_filepath = "Data/plate1.final_data.Gompertz_with_bump.tsv"
    QV_df = pd.read_csv(QV_filepath, sep="\t", usecols=["WellID", "Sigmoid_curve", "Condition", "Gene", "Theoretical_N0_log"] + cycles_col)
    QV_df = QV_df.loc[(QV_df.Sigmoid_curve) & QV_df.Condition.isin(["PS" + str(x) for x in range(1, 13)]), :].copy()
    # Filtering out gene/condition without the 6 replicates showing a sigmoid curve
    QV_df = QV_df.groupby(["Gene", "Condition"]).filter(lambda row: row["WellID"].count() == 6).copy()
    # Loading research dataset
    research_filepath = "Data/Plate2.final_data.Gompertz_with_bump.REF.tsv"
    research_df = pd.read_csv(research_filepath, sep="\t", usecols=["WellID", "Sigmoid_curve", "Pre_amplification", "Condition", "Gene", "Theoretical_N0_log"] + cycles_col)
    research_df = research_df.loc[(research_df.Sigmoid_curve) & (research_df.Pre_amplification == 0), :].copy()
    # Filtering out gene/condition without the 6 replicates showing a sigmoid curve
    research_df = research_df.groupby(["Gene", "Condition"]).filter(lambda row: row["WellID"].count() == 6).copy()
    # Setting the training/testing plate and validation one
    if options.config1:
        ttp = research_df
        excluded_gene_condition_ttp = [("Notch1a", "PC"), ("Notch1a", "PS"), ("BMP4", "PC"), ("Crabp2a", "PC")]
        vp = QV_df
        excluded_gene_condition_vp = [("Bactine2", "PS")]
    else:  # Config2
        ttp = QV_df
        excluded_gene_condition_ttp = [()]
        vp = research_df
        excluded_gene_condition_vp = [("Notch1a", "PC"), ("Notch1a", "PS"), ("BMP4", "PC"), ("Crabp2a", "PC"), ("Bactine2", "PS"), ("Bactine2", "PC")]  # To validate on the research plate, we discard BActine2 to be fair with other metohds (used as reference for Cq method)

    # Loading the train/test dataset
    wellID_ttp, (gene_train, condition_train, x_train, y_train), (gene_test, condition_test, x_test, y_test), my, sy = read_file(ttp, test_dataset=test_gene_condition, exclude=excluded_gene_condition_ttp)
    # Creating the deep neural network
    model = Sequential()  # Initialization
    model.add(InputLayer(input_shape=(40, 1)))  # Providing the shape of our data
    # Chaining layer transformations on top of the input
    model.add(Conv1D(5, kernel_size=5, strides=1, activation="sigmoid"))  # Convolutional layer with "kernel_size" for the size of the kernel to slide, "strides" for the steps size of the sliding window
    model.add(MaxPooling1D(pool_size=2, strides=2))  # Smoothing data from values around
    model.add(Flatten())
    model.add(Dense(3, activation="sigmoid"))  # Summarizing data to a vector of 3
    model.add(Dense(1))  # Getting result as one value

    # Printing the summary of how the data will get transformed at each stage of the model
    model.summary()

    # Compiling the model
    model.compile(loss=keras.losses.mse,  # Mean Squared Error
                  optimizer=keras.optimizers.Adam(),  # A stochastic gradient descent method based on adaptive estimation of first-order and second-order moments
                  metrics=["mae"])  # Mean Absolute Error
    # Fitting the model to the data
    nb_epoch = 2500
    history = model.fit(x_train, y_train,
              batch_size=128,  # Takes curves by batch of 128
              epochs=nb_epoch,
              verbose=2,
              validation_data=(x_test, y_test))

    # Plotting loss and metric
    for k in history.history.keys():
        plt.plot(range(nb_epoch), history.history[k], label=k)
        plt.legend()
    plt.show()

    # Displaying the correlation between theoretical and predicted for the train set
    display_plot(model, x_train, y_train, my, sy, gene_train)
    # Displaying the correlation between theoretical and predicted for the test set
    display_plot(model, x_test, y_test, my, sy, gene_test)

    # Validating the model on the other plate
    wellID_vp, (gene_val, condition_val, x_val, y_val), my, sy = read_file(vp, validation=True, exclude=excluded_gene_condition_vp)
    # Displaying the correlation between theoretical and predicted for the validation set
    display_plot(model, x_val, y_val, my, sy, gene_val)
    # Reporting the predicted values
    if options.config1:
        np.savetxt("DL_predictions.plate1.tsv", np.column_stack((wellID_vp, back_convert(model.predict(x_val).T[0], my, sy))), delimiter="\t", fmt="%s")
    else:  # Config2
        np.savetxt("DL_predictions.plate2.tsv", np.column_stack((wellID_vp, back_convert(model.predict(x_val).T[0], my, sy))), delimiter="\t", fmt="%s")
