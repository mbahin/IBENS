#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Common functions for qPCR scripts. """

import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Setting global variables
cycles_col = ["F" + str(i).zfill(2) for i in range(40)]
genes = ["BMP4", "Bactine2", "Crabp2a", "Cyp26a1", "DeltaC", "DeltaD", "FGF8a", "Her7", "Hoxa1a", "Hoxa5a", "Hoxb1b", "Hoxd4a", "Notch1a", "RPL13a1", "RPL13a2", "Rarab"]


def decorticate_condition(condition):
    ''' To decompose the condition into pre-amplification cycle number, PS or PC and dilution factor. '''
    match = re.match(r'((1[048])C_)?(P[CS])([1-9][0-2]?)', condition)
    # For the plate2, there will always be a match but not for plate1
    if match:
        # Checking if there was pre-amplification cycles
        if not match.group(2):
            pre_amp = 0
        else:
            pre_amp = int(match.group(2))
        return pre_amp, match.group(3), int(match.group(4))
    else:
        return np.nan, np.nan, np.nan


def compute_N0(row, molecules):
    # If the well is water or if it was pre-amplified, then there is no theoretical N0
    if (row.Dilution == 12) or (row.Pre_amplification != 0):  # Warning, in plate2, the 12th dilution are water samples, not in plate2 but always flat curves so it's ok
        return np.nan
    else:
        # Applying the dilution factor to the known starting molecule quantity
        return molecules.loc[str(row.S_or_C) + "1", row.Gene] * np.power(10, 1 - float(row.Dilution)) * 0.0027
        """
        n = N0_df.loc[(N0_df.Gene == row.Gene) & (N0_df.Condition == str(row.S_or_C) + "1"), "Molecules_per_µl"].iloc[0] * np.power(10, 1 - float(row.Dilution))
        # Applying different formula if there was pre-amplification or not
        if row.Pre_amplification == 0:
            return n * 0.0027
        else:
            return n * 0.000135 * np.power(E, float(row.Pre_amplification))
        """


def gather_data(filepath, sample_size, water_samples, output_filepath, N0_file):
    print("Loading qPCR data from file " + filepath + "...")
    # Loading intermediate data
    df = pd.read_csv(filepath, sep="\t", nrows=sample_size)

    # Defining metadata (Warning: in plate1, Preamp and Dilution are float, probably because of NaNs)
    df[["Pre_amplification", "S_or_C", "Dilution"]] = df.apply(lambda row: decorticate_condition(row.Condition), axis=1, result_type="expand")

    # Discriminating and discarding flat curves from the dataset
    df = discriminateFlatCurves_pd(df, sample_size, water_samples, output_filepath)
    # Determining the rising cycle
    df["Rising_cycle"] = df.apply(lambda row: 0 if not row["Sigmoid_curve"] else area_under_chord_pd(row[cycles_col], 11)[1], axis=1)

    # Loading the theoretical N0 for PS1 (and PC1 if plate2)
    molecules = pd.read_csv(N0_file, sep="\t", header=0, index_col=0)
    # Computing the theoretical N0
    df["Theoretical_N0"] = df.apply(lambda row: compute_N0(row, molecules), axis=1)
    df["Theoretical_N0_log"] = np.log10(df["Theoretical_N0"])

    return df


def read_fluo_raw_values(filepath, nb_rows_to_skip):
    ''' Reads a block of 9.216 x 40 raw fluo values (either FAM, ROX, background ROX or background FAM). '''
    # Loading the 9.216 x 40 values
    df = pd.read_csv(filepath, sep=",", usecols=range(41), skiprows=nb_rows_to_skip, nrows=9216, header=None)
    # Renaming the columns
    df.columns = ["WellID"] + list(range(40))
    # Defining the well ID as index
    df.set_index("WellID", inplace=True)
    return df


def extract_data_from_Fluidigm_file(filepath, output_filepath):
    ''' Extracts the raw data from Fluidigm file and compute the fluo values (from FAM, ROX and their background) to create the HiFit input format. '''
    print("Loading data from Fluidigm raw file...")
    # Loading well infos
    infos = pd.read_csv(filepath, sep=",", names=["WellID", "Condition", "Gene", "Fluidigm_Ct"], usecols=[0, 1, 4, 6], skiprows=12, nrows=9216)
    # Creating the sample and assay columns
    infos["Sample"] = infos.apply(lambda row: int(row.WellID[1:3]), axis=1)
    infos["Assay"] = infos.apply(lambda row: int(row.WellID[5:]), axis=1)
    # Loading the 4 raw fluo values
    ROX = read_fluo_raw_values(filepath, 9236)
    FAM = read_fluo_raw_values(filepath, 18455)
    FAM.to_csv("/data/biocomp/bahin/FAM.tsv", sep="\t")
    Bkgd_ROX = read_fluo_raw_values(filepath, 27674)
    Bkgd_FAM = read_fluo_raw_values(filepath, 36893)
    # Computing the fluo values, renaming the columns and setting the well ID as a column
    fluo = (FAM - Bkgd_FAM) / (ROX - Bkgd_ROX)
    fluo.columns = ["F" + str(x).zfill(2) for x in range(40)]
    fluo.reset_index(inplace=True)
    # Merging well infos and fluo values
    df = infos.merge(fluo, on="WellID")
    # Storing the output
    #df.to_csv(output_filepath, index=False, sep="\t")
    return df


def area_under_chord_pd(data, window_size=11):  ## Problem of minus value (because of the odd window size but add +1 value left or right??
    ''' Function computing the area under the curve for the data and according to a defined window size. Its returns the computed areas under the chord, the data point with the maximum area and its value. '''
    # Checking window_size oddness
    if window_size % 2 == 0:
        sys.exit("Error in 'area_under_chord' function, an even number was given as window size. Aborting.")
    # Initializing the areas array with zeros for half the size of the window (so that the area value is central for the window in the array)
    areas = np.array([np.zeros(int(window_size / 2))])
    # Transforming the fluo values list onto a Series (called from discriminateFlatCurves it's already a Series)
    data = pd.Series(data)
    # Looping over the windows
    for window_start in range(0, len(data) - window_size):
        # Getting the linear regression coordinates for the considered window (not using stats.linregress because for 2 points, r = 1.0 and it raises a numpy warning)
        slope = (data.iloc[window_start + window_size] - data.iloc[window_start]) / window_size
        intercept = data.iloc[window_start] - (slope * window_start)
        # Computing the area under the chord (sum of the differences between the chord and the curve over the window length)
        areas = np.append(areas, np.sum([(((window_start + p) * slope) + intercept) - data.iloc[window_start + p] for p in range(window_size)]))
    # The result is the central point of the window with the maximum area under the chord or 1
    areas = np.hstack((areas, np.zeros(int(window_size / 2))))
    return areas, np.argmax(areas), areas[np.argmax(areas)]


def gaussian(x, mu, sig):
    ''' Return data to draw a Gaussian curve. '''
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def discriminateFlatCurves_pd(df, sample_size, water_samples, output_filepath):
    ''' Function to discard the flat curves. '''
    # Computing curve amplitudes (max - min difference)
    df["Amplitude"] = df.apply(lambda row: row[cycles_col].max() - row[cycles_col].min(), axis=1)
    # Sorting the amplitudes
    amplitudes_order = df["Amplitude"].sort_values()

    # Plotting the measures and the area under the chord
    results_AUChord = area_under_chord_pd(amplitudes_order, 11)
    derv_measures = amplitudes_order.diff()  # Computing the derivative (difference between a value and the next one)
    derv_measures.iat[0] = 0  # Fixing the first NaN value to 0
    g = gaussian(np.linspace(0, 1, sample_size), 0.5, 0.25)  # To create gaussian distribution
    """
    fig = plt.figure(figsize=(8, 5))  # To create plt object
    ax1 = plt.gca()  # Setting the first axis
    ax2 = ax1.twinx()  # Building the second axis
    ax1.plot(range(amplitudes_order.size), amplitudes_order, color="#1f77b4", label="Sorted amplitudes")  # Plotting the area under the chord distribution
    ax1.plot(g, color="#c7c7c7", label="Gaussian distribution")  # Plotting the Gaussian distribution
    ax2.plot(results_AUChord[0])  # Plotting the area under the chord distribution
    ax2.plot(range(amplitudes_order.size), derv_measures, color="#d62728", label="Dervative")  # Derivatives
    ax2.plot(range(amplitudes_order.size), derv_measures * g, color="#2ca02c", label="Revised derivative")  # Plotting the "transformed" area under the chord distribution
    ax1.legend(loc="upper left", shadow=True)
    ax2.legend(loc="upper center", shadow=True)
    ax1.set_ylim(0, 1.2)
    ax2.set_ylim([0, 0.25])
    ax2.set_axis_off()
    plt.savefig(output_filepath + "flat_curves.png")
    print("Threshold curve:", df.iloc[(derv_measures * g).idxmax()])
    """
    # Getting the threshold between flat and sigmoid curves and assessing curves status
    amplitude_threshold = df.iloc[(derv_measures * g).idxmax()].Amplitude
    df["Sigmoid_curve"] = df.apply(lambda row: row.Amplitude > amplitude_threshold, axis=1)

    # Setting water samples as flat
    df.loc[df.Sample.isin(water_samples), "Sigmoid_curve"] = False

    # Printing stats
    discarded_wells = len(df.index) - df["Sigmoid_curve"].sum()  # Counts the number of rows with "Sigmoid_curve" True
    print("Out of the " + str(sample_size) + " considered curves, " + str(discarded_wells) + " are found to be flat (amplitude < " + str(amplitude_threshold) + ") or from water samples. They will be discarded for the further steps (" + str(sample_size - discarded_wells) + " remaining curves).")

    # Creating a flat/sigmoid curves heatmap
    """
    heatmap = pd.pivot_table(df, values="Sigmoid_curve", index="Sample", columns="Assay")
    fig = plt.figure(figsize=(10, 10))
    ax = sns.heatmap(heatmap, xticklabels=6, yticklabels=12, cmap="RdYlBu", linewidths=0.05, cbar=False, square=True)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.title("Plate sigmoid/flat curves")
    plt.savefig(output_filepath + "sigmoid_flat_curves_heatmap.png")
    """
    """
    # Reporting the sample IDs for N random flat curves, the Nth almost sigmoid flat curves, the Nth almost flat sigmoid curves and N random sigmoid curves
    n = 25
    sample_IDs = data_infos_global[:, 0:1]
    with open("flat.random.txt", "w") as output_file:
        for _ in range(n):
            #output_file.write(np.random.choice(sample_IDs[measures_order][0:flat_curves_threshold])[0] + "\n")
            output_file.write(np.random.choice(sample_IDs[measures_order][0:flat_curves_threshold].flatten()) + "\n")
    with open("flat.edge.txt", "w") as output_file:
        for i in range(n):
            output_file.write(sample_IDs[measures_order][(flat_curves_threshold - i - 1):flat_curves_threshold][0][0] + "\n")
    with open("sigmoid.random.txt", "w") as output_file:
        for _ in range(n):
            #output_file.write(random.choice(sample_IDs[measures_order][flat_curves_threshold:])[0] + "\n")
            output_file.write(np.random.choice(sample_IDs[measures_order][flat_curves_threshold:].flatten()) + "\n")
    with open("sigmoid.edge.txt", "w") as output_file:
        for s in sample_IDs[measures_order][flat_curves_threshold:(flat_curves_threshold + n):]:
            output_file.write(s[0] + "\n")
        #for i in range(n):
            #output_file.write(sample_IDs[measures_order][flat_curves_threshold:(flat_curves_threshold + i + 1)][0][0] + "\n")
    """

    return df


def correct_bump(row, bump1, bump2, background=0):
    ''' Returns the fluo values from which the modelled bump (Gompertz_with_bump equation) + possibly the background is subtracted. '''
    corr = [bump1 - bump1 * np.exp(x * bump2) + background for x in range(40)]
    return row[cycles_col] - corr

"""
def create_N0_df(N0_file, E=2):
    ''' [OLD: we can't really use preamplified data] Creates the theoretical N0 DataFrame (depends on the genes efficiency assumed). '''

    # Loading PS1 and PC1 values for each gene
    molecules = pd.read_csv(N0_file, sep="\t", header=0, index_col=0)
    # Defining the samples
    pre_amp = ["", "10C_", "14C_", "18C_"]  # With the not preamplified condition because we need it as a base to compute the pre-amplification
    samples = [i + j + str(k) for i in pre_amp for j in ["PS", "PC"] for k in range(1, 13)]
    # Creating the DataFrame skeleton
    N0_df = pd.DataFrame(index=samples, columns=genes, dtype=float)
    N0_df.index.names = ["Condition"]
    # Copying the input PS1 and PC1 molecules count into the DataFrame
    N0_df.loc["PS1"] = molecules.loc["PS1"]
    N0_df.loc["PC1"] = molecules.loc["PC1"]
    # Melting the DataFrame to have one row per condition and gene
    N0_df = N0_df.reset_index().melt(id_vars="Condition", var_name="Gene", value_name="Molecules_per_µl")
    N0_df[["Pre_amplification", "S_or_C", "Dilution"]] = N0_df.apply(lambda row: decorticate_condition(row.Condition), axis=1, result_type="expand")
    N0_df["Theoretical_N0"] = N0_df.apply(lambda row: compute_N0(row, N0_df, E), axis=1)
    N0_df["Theoretical_N0_log"] = np.log10(N0_df["Theoretical_N0"])
    return N0_df
"""
