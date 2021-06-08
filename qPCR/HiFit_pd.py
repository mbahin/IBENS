#!/usr/bin/python
# -*- coding: utf-8 -*-

# This Python script intended to do absolute quantification of microfluidics qPCR array using sigmoidal model.
# (We want to find the initial DNA amount of a qPCR sample using a single curve.)

import sys
import argparse
import re
import string
from random import random
import numpy as np
from scipy.optimize import minimize
from sklearn.cluster import KMeans
import pandas as pd
from dask import dataframe as dd
from dask.multiprocessing import get
from dask.diagnostics import ProgressBar
# Loading modules
import qPCR_functions

###################################################################################################
# Functions

def logistic4p(x, minimum, maximum, F0, E):
    ''' Equation of a four parameters sigmoid curve (Verhulst model)
    Input: "x" a integer or an np.array of integer(the index of a point of a qPCR curve), minimum/maximum/F0/E a float
    Output: f(x) a float or an array of float, the image of x by the equation
    '''
    try:
        return minimum + (maximum - minimum) / (1 + ((((maximum - minimum) / F0) - 1) * np.power(E, -x)))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([0], len(x))


def logistic4p_with_bump(x, minimum, maximum, F0, E, bump1, bump2):
    ''' Equation of a four parameters sigmoid curve (Verhulst model)
    Input: "x" a integer or an np.array of integer(the index of a point of a qPCR curve), minimum/maximum/F0/E a float
    Output: f(x) a float or an array of float, the image of x by the equation
    '''
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + (maximum - minimum) / (1 + ((((maximum - minimum) / F0) - 1) * np.power(E, -x)))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def logistic5p(x, minimum, maximum, IP, slope, f):
    try:
        return minimum + (maximum - minimum) / (np.power(1 + np.exp(slope * (x - IP)), f))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def logistic5p_with_bump(x, minimum, maximum, IP, slope, f, bump1, bump2):
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + (maximum - minimum) / (np.power(1 + np.exp(slope * (x - IP)), f))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def logistic5p_with_bump_log(x, minimum, maximum, IP, slope, f, a, b):  ## Bump not well modelled (not like others)
    try:
        return a * np.exp(b * x) + a + minimum + (maximum - minimum) / (np.power(1 + np.exp(slope * (np.log(x) - np.log(IP))), f))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used)
        return np.repeat([np.nan], len(x))


def logistic6p(x, minimum, maximum, F0, E, alpha, beta):
    ''' Equation of a six parameters sigmoid curve
    Input: "x" a integer or an np.array of integer(the index of a point of a qPCR curve), minimum/maximum/F0/E/alpha/beta a float
    Output: f(x) a float or an array of float, the image of x by the equation
    '''
    try:
        return minimum + (maximum - minimum) / (1 + ((((maximum - minimum) / F0) - 1) * np.power(E, -x)) + alpha * np.power((F0 * np.power(E, x)), -beta))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def logistic6p_with_bump(x, minimum, maximum, F0, E, alpha, beta, bump1, bump2):
    ''' Equation of a six parameters sigmoid curve
    Input: "x" a integer or an np.array of integer(the index of a point of a qPCR curve), minimum/maximum/F0/E/alpha/beta a float
    Output: f(x) a float or an array of float, the image of x by the equation
    '''
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + (maximum - minimum) / (1 + ((((maximum - minimum) / F0) - 1) * np.power(E, -x)) + alpha * np.power((F0 * np.power(E, x)), -beta))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def Gompertz(x, a, b, c, d):
    try:
        return a + (b - a) * np.exp(-1 * np.exp(-d * (x - c)))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def Gompertz_with_bump(x, minimum, maximum, c, d, bump1, bump2):
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + (maximum - minimum) * np.exp(-1 * np.exp(-d * (x - c)))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def MAK2(x, minimum, D0, k):
    if x == 0:
        return D0
    else:
        return (x - 1) + k * np.log(1 + ((x -1) / k))


def Boltzmann_LRE_with_bump(x, minimum, maximum, IP, slope, bump1, bump2):
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + maximum / (1 + np.exp(-1 * ((x - IP) / slope)))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))


def scf(x, Fb, Fmax, c, b):
    try:
        return (Fmax / (1 + np.exp(((-1) / b) * (x - c)))) + Fb
    except FloatingPointError:
        return np.repeat([0], len(x))


def scf_with_bump(x, minimum, maximum, c, b, bump1, bump2):
    try:
        return bump1 - bump1 * np.exp(x * bump2) + (maximum / (1 + np.exp(((-1) / b) * (x - c)))) + minimum
    except FloatingPointError:
        return np.repeat([np.nan], len(x))


def Richards(x, Fb, Fmax, c, b, d):
    try:
        return Fb + (Fmax / np.power(1 + np.exp((c - x) / b), d))
    except FloatingPointError:
        return np.repeat([0], len(x))


def Richards_with_bump(x, minimum, maximum, c, b, d, bump1, bump2):
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + (maximum / np.power(1 + np.exp((c - x) / b), d))
    except FloatingPointError:
        return np.repeat([np.nan], len(x))


def Hill_with_bump(x, minimum, maximum, c, d, bump1, bump2):
    try:
        return bump1 - bump1 * np.exp(x * bump2) + minimum + (maximum * np.power(x, c)) / (np.power(d, c) + np.power(x, c))
    except FloatingPointError:
        return np.repeat([np.nan], len(x))


def residualSumSquare(parameter_rss, row, equation):
    ''' Computes the residual sum square.
    Input: a np.array "parameter_rss" containing parameters fitting the sigmoid curve from "data_rss" points, a list cycles on which to compute the RSS and en equation.
    Output: the residual sum square (square difference between the actual curve and the fitted one).
    '''
    if parameter_rss[0] > parameter_rss[1]:  # The minimize function with "SLSQP" method sometimes fails to consider the constraints (https://github.com/scipy/scipy/pull/8986)
        return np.nan
    """ Faster version with fluo values in a list in one column
    fitting_start_cycle = max(int(row.Fluidigm_Ct) - 5, 0)
    return np.sum((row.Fluo[fitting_start_cycle:] - equation(np.array(range(fitting_start_cycle, 40)), *parameter_rss))**2)
    """
    try:
        return np.sum((np.array(row[cycles_col]) - equation(np.arange(40), *parameter_rss))**2)
    except (OverflowError, FloatingPointError):
        return np.nan

###################################################################################################
# HiFit model

def onefit(row, equation):
    ''' Function that tries to optimize the fitting of a curve by a model.
    Inputs: the sample ID, the fluo values of the sample, the IPS to use and the equation model.
    Outputs: the sample ID, the final parameter set, the cycles on which the fitting was optimized and the final RSS.
    '''

    # Defining a constraint for the minimum and maximum parameters to be different
    cons = {"type": "ineq", "fun": lambda x: x[1] - x[0]}

    # Setting the parameters bounds for the fitting equation
    params = ["P_" + letter + "_" + equation.__name__ for letter in list(string.ascii_lowercase)[0:equation.__code__.co_argcount - 1]]
    minimized_IPS = minimize(residualSumSquare, x0=row[params], args=(row, equation), bounds=bounds[equation.__name__], constraints=cons, method="SLSQP")
    if np.isnan(minimized_IPS.x[0]):  # In small tests, usually, all solutions are bad so we return a fake params and a RSS of 1 not to crash
        return np.repeat([1], len(params) + 1)
    else:
        RSS = residualSumSquare(minimized_IPS.x, row, equation)
        return np.append(minimized_IPS.x, RSS)


def fit(data, equation, nb_init_param_sets, nb_processors, pinit=None):
    ''' Function that multiprocessingly fits curves according to initial parameters sets. '''
    # Creating a random sets of initial parameters
    if not isinstance(pinit, pd.DataFrame):
        pinit = pd.DataFrame()
        for letter, bound in zip(list(string.ascii_lowercase)[0:equation.__code__.co_argcount - 1], bounds[equation.__name__]):
            pinit["P_" + letter + "_" + equation.__name__] = pd.Series(np.random.uniform(bound[0], bound[1], nb_init_param_sets))
        # Swapping the min and max parameters if the second (end level) is lesser than the first (start level)
        pinit.loc[pinit["P_a_" + equation.__name__] > pinit["P_b_" + equation.__name__], ["P_a_" + equation.__name__, "P_b_" + equation.__name__]] = pinit.loc[pinit["P_a_" + equation.__name__] > pinit["P_b_" + equation.__name__], ["P_b_" + equation.__name__, "P_a_" + equation.__name__]].values

    # Creating a n*m DataFrame (with the n curves from the subset * the m random parameter sets)
    data["key"] = 0
    pinit["key"] = 0
    initialization_matrix = data.merge(pinit, how="outer")

    # For each row of the new DataFrame, applying a fitting to the fluo data from the initial parameter sets
    final_params = ["P_" + letter + "_final_" + equation.__name__ for letter in list(string.ascii_lowercase)[0:equation.__code__.co_argcount - 1]] + ["RSS_" + equation.__name__]
    initialization_matrix[final_params] = dd.from_pandas(initialization_matrix, npartitions=nb_processors).map_partitions(lambda df: df.apply(lambda row: onefit(row, equation), axis=1, result_type="expand")).compute(scheduler=get)

    # For each curve, keep the IPS minimizing the RSS
    initialization_matrix.dropna(inplace=True)  # For the little test, if one well has only NaN value
    optimized_params = initialization_matrix.loc[initialization_matrix.groupby("WellID")["RSS_" + equation.__name__].idxmin()] # If there are several solutions with the same RSS, a random one among them is taken
    return optimized_params


def initialFitting(init_size, equation, nb_init_param_sets, stable_test, nb_param_sets_clusters, nb_processors):
    ''' Function that finds a set of optimized initial parameters sets for sigmoid curves within the dataset. '''
    print("Trying " + str(nb_init_param_sets) + " random initial parameter sets (IPS) on a subset of " + str(init_size) + " curves...")

    # Randomly picking the subset of curves that will be used for the initialization
    if stable_test:
        curves_subset = data.sample(init_size, random_state=45)
    else:
        curves_subset = data.sample(init_size)

    # Running the fitting model on the subsets of curves
    curves_subset = fit(curves_subset, equation, nb_init_param_sets, nb_processors)
    final_params = ["P_" + letter + "_final_" + equation.__name__ for letter in list(string.ascii_lowercase)[0:equation.__code__.co_argcount - 1]]
    optimized_params = curves_subset[final_params]

    # Clustering the optimized parameters
    print("Computing the " + str(nb_param_sets_clusters) + " bests final parameter sets by clustering...")
    # Rescaling parameters to comparable values
    optimized_params_low = optimized_params.quantile(0.05)
    optimized_params_high = optimized_params.quantile(0.95)
    optimized_params_rescaled = (optimized_params - optimized_params_low.values) / (optimized_params_high.values - optimized_params_low.values)
    # Clustering the initial parameter sets
    km = KMeans(nb_param_sets_clusters)
    try:
        clusters = km.fit(optimized_params_rescaled)
    except ValueError:
        print(optimized_params_rescaled)
    kcenters = pd.DataFrame(clusters.cluster_centers_)
    kcenters = pd.DataFrame(clusters.cluster_centers_)
    # kmeans_labels = a.labels_  # To visualize the clusters
    # Rescaling the k-centers to to initial range of values
    optimized_params_clustered = (kcenters * (optimized_params_high.values - optimized_params_low.values)) + optimized_params_low.values
    # Renaming the parameter columns
    optimized_params_clustered.columns = ["P_" + letter + "_" + equation.__name__ for letter in list(string.ascii_lowercase)[0:equation.__code__.co_argcount - 1]]

    return optimized_params_clustered


def HiFitModel(data, init_size, equation, nb_init_param_sets, stable_test, nb_param_sets_clusters, nb_processors, IPS_file, output_filepath):

    print("Running HiFit model with " + str(equation.__name__) + " equation...")
    # Checking whether initial parameter sets is provided and with the appropriate parameters
    if (IPS_file == "") or ("P_a_" + equation.__name__ not in pd.read_csv(IPS_file, sep="\t", index_col=0, header=0).columns):
        print("Creating exhaustive initial parameter sets")
        initial_parameters = initialFitting(init_size, equation, nb_init_param_sets, stable_test, nb_param_sets_clusters, nb_processors)
        if IPS_file != "":
            # Updating the provided initial parameter sets file
            pd.concat([pd.read_csv(IPS_file, sep="\t", index_col=0, header=0), initial_parameters], axis=1).to_csv(IPS_file, sep="\t", index=False)
        else:
            # Creating a new initial parameter sets file
            initial_parameters.to_csv(output_filepath + "/IPS.tsv", sep="\t", index=False)
            # Setting the new IPS file path for the next equation if there is one
            IPS_file = output_filepath + "/IPS.tsv"
    else:
        print("Using previously created initial parameter sets")
        initial_parameters = pd.read_csv(IPS_file, sep="\t", index_col=0, header=0)

    # Fitting all data quickly thanks to the previous step
    print("Smart fitting " + equation.__name__ + "...")
    data = fit(data, equation, nb_init_param_sets, nb_processors, initial_parameters)

    # Creating the modelled fluo values (with min and max to 0 and 1)
    #mod_cycles_col = ["F" + str(i).zfill(2) + "_mod_" + equation.__name__ for i in range(40)]
    #data[mod_cycles_col] = data.apply(lambda row: modelled_fluo(row, equation), axis=1, result_type="expand")

    # Creating a DataFrame with one modelled fluo value per row
    #mod_fluo = pd.melt(data, id_vars=["Gene", "WellID", "Condition", "Pre_amplification", "S_or_C", "Dilution"], var_name="Cycle", value_vars=mod_cycles_col, value_name="Modelled_fluo_value")
    # Removing the initial "F" and final "_mod" of modelled cycles column names
    #mod_fluo["Cycle"] = mod_fluo.apply(lambda row: int(row["Cycle"][1:-4]), axis=1)

    data.to_csv(output_filepath + "/final_data." + equation.__name__ + ".tsv", sep="\t", index=False)
    return data, IPS_file

###################################################################################################
# Main
if __name__ == "__main__":
    # All numpy exceptions will be raised
    np.seterr(all="raise")

    # Setting Dask progress bar
    ProgressBar().register()

    # Setting the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--Fluidigm_file", dest="Fluidigm_file", help="Fluidigm raw output file")
    parser.add_argument("-f2", "--intermediate_file", dest="HiFit_intermediate_file", help="HiFit intermediate file (data read from Fluidigm raw file but not yet gathered)")
    parser.add_argument("-d", "--data", type=str, dest="data_file", default="Data/plate2.tsv", help="Path to the input data file to analyse.", metavar="FILE")
    parser.add_argument("-b", "--bump", action="store_true", dest="bump", help="Whether the bump hasn't yet been removed from the data.")
    parser.add_argument("-r", "--research_plate", action="store_true", dest="research_plate", help="Whether this is THE research plate.")
    parser.add_argument("-m", "--model", dest="model", nargs="+", required=True, help="Model to use for the fitting.")
    parser.add_argument("-s", "--sample-size", type=int, dest="sample_size", default=9216, help="Number of rows of the qPCR plate to analyze.(default= 9216 (96*96))")
    parser.add_argument("-i", "--init-size", type=int, dest="init_size", default=1000, help="Number of random curves to create the initial parameter sets.(default 1k)")
    parser.add_argument("-n", "--nb-init-param-sets", type=int, dest="nb_init_param_sets", default=100, help="Number of random initial parameter sets to create the initial parameter sets.(default=100)")
    parser.add_argument("-c", "--nb-param_sets-clusters", type=int, dest="nb_param_sets_clusters", default=20, help="Number of clusters of initial parameter sets to create.(default=20)")
    parser.add_argument("-t", "--stable-test", action="store_true", dest="stable_test", help="Flag to specify that we use the stable test.")
    parser.add_argument("-l", "--little-test", action="store_true", dest="little_test", help="Flag to specify that we use the little test.")
    parser.add_argument("-p", "--nb_processors", type=int, dest="nb_processors", default=70, help="Number of processors used by the script.(default=4)")
    parser.add_argument("--IPS_file", dest="IPS_file", default="", help="Providing an IPS file to ditch the first phase")
    parser.add_argument("-o", "--output", dest="output", help="Output filepath")
    options = parser.parse_args()

    # Setting fixed parameters for the stable test
    if options.stable_test:
        options.sample_size = 960
        options.init_size = 100
        options.nb_init_param_sets = 10
        options.nb_param_sets_clusters = 3
        #options.label = "stable"
        np.random.seed(45)  # Fixing an arbitrary seed for reproducibility
    elif options.little_test:
        options.sample_size = 384
        options.init_size = 10
        options.nb_init_param_sets = 5
        options.nb_param_sets_clusters = 3
    else:
        # Ensuring that the sample size is greater than or equal to 96 and is a multiple of 96
        if options.sample_size < 96:
            options.sample_size = 96
        else:
            options.sample_size = int(options.sample_size / 96) * 96
        # Ensuring that the init size is lesser than or equal to the sample size
        if options.init_size > options.sample_size:
            options.init_size = options.sample_size

    # Setting fluo values to process
    if options.bump:
        # If the bump hasn't yet been removed
        cycles_col = ["F" + str(i).zfill(2) for i in range(40)]
    else:
        # For the data with the bump removed
        cycles_col = ["F" + str(i).zfill(2) + "_bump_corr" for i in range(40)]  ## TO MOVE TO qPCR_functions.py if validated with AG

    # Setting the dict to match the model argument (string) with the actual function
    models = {"logistic4p": logistic4p,
              "logistic4p_with_bump": logistic4p_with_bump,
              "logistic5p": logistic5p,
              "logistic5p_with_bump": logistic5p_with_bump,
              "logistic5p_with_bump_log": logistic5p_with_bump_log,
              "logistic6p": logistic6p,
              "logistic6p_with_bump": logistic6p_with_bump,
              "Gompertz": Gompertz,
              "Gompertz_with_bump": Gompertz_with_bump,
              "Boltzmann_LRE_with_bump": Boltzmann_LRE_with_bump,
              "scf": scf,
              "scf_with_bump": scf_with_bump,
              "Richards": Richards,
              "Richards_with_bump": Richards_with_bump,
              "Hill_with_bump": Hill_with_bump,
              }
    # Setting variables
    bounds = {"logistic4p": ((0, 1), (0, 2), (1 * 10 ** (-8), 1.), (0.5, 2)),
              "logistic4p_with_bump": ((0, 0.5), (0, 2), (10 ** (-8), 1.), (1, 10), (0, 0.5), (-0.5, 0)),
              "logistic5p": ((0, 0.5), (0, 2), (0, 40), (-1, 0), (0, 50)),
              "logistic5p_with_bump": ((0, 0.5), (0, 2), (0, 40), (-1, 0), (0, 50), (0, 0.5), (-0.5, 0)),
              "logistic5p_with_bump_log": ((0, 1), (0, 2), (0, 40), (-30, -5), (0, 30), (-3, 1), (-10, 1)),
              "logistic6p": ((0, 1), (0, 2), (1 * 10 ** (-8), 10), (0, 5), (-50, 50), (-100, 100)),
              "logistic6p_with_bump": ((0, 0.5), (0, 2), (1 * 10 ** (-8), 1.), (0, 3), (-30, 10), (-100, 100), (0, 0.5), (-0.5, 0)),
              "Gompertz": ((0, 0.5), (0, 2), (0, 40), (0, 1)),
              "Gompertz_with_bump": ((0, 0.5), (0, 2), (0, 40), (0, 1), (0, 0.5), (-0.5, 0)),
              "Boltzmann_LRE_with_bump": ((0, 0.5), (0, 2), (0, 40), (0, 10), (0, 0.5), (-0.5, 0)),
              "scf": ((0, 1), (0, 2), (0, 40), (0, 5)),
              "scf_with_bump": ((0, 0.5), (0, 2), (0, 40), (0, 5), (0, 0.5), (-0.5, 0)),
              "Richards": ((0, 0.5), (0, 2), (0, 40), (0, 40), (-1, 30)),
              "Richards_with_bump": ((0, 0.5), (0, 2), (0, 40), (0, 40), (-1, 30), (0, 0.5), (-0.5, 0)),
              "Hill_with_bump": ((0, 0.5), (0, 2), (0, 30), (0, 40), (0, 0.5), (-0.5, 0))}
    # Setting water samples
    if options.research_plate:
        water_samples = [12, 24, 36, 48, 60, 72, 84, 96]
    else:  # Plate1
        water_samples = [1, 13, 25, 37, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84]

    #################################################

    # Extracting the data from Fluidigm raw file if it was provided
    if options.Fluidigm_file:
        qPCR_functions.extract_data_from_Fluidigm_file(options.Fluidigm_file, options.HiFit_intermediate_file)

    # Gathering the data if the intermediate file is provided
    if options.HiFit_intermediate_file:
        # Setting N0 file
        if options.research_plate:
            N0_file = "Data/plate2.PS1_PC1_N0.tsv"
        else:  # Plate1
            N0_file = "Data/plate1.PS1_N0.tsv"
        # Loading fluo and metadata
        data = qPCR_functions.gather_data(options.HiFit_intermediate_file, options.sample_size, water_samples, options.output, N0_file)
        data.to_csv(options.data_file, sep="\t", index=False)
    else:  # Reading the gathered data directly
        data = pd.read_csv(options.data_file, sep="\t", nrows=options.sample_size)
        data = data.sample(options.sample_size)

    # Discarding the flat curves
    data = data.loc[data.Sigmoid_curve]

    # Ensuring that the initial size is lesser than or equal to the number of sigmoid curves from the sample size
    init_size = min(len(data.index), options.init_size)

    # Running HiFit for the provided models
    IPS_file = options.IPS_file
    for model in options.model:
        data, IPS_file = HiFitModel(data, init_size, models[model], options.nb_init_param_sets, options.stable_test, options.nb_param_sets_clusters, options.nb_processors, IPS_file, options.output)
