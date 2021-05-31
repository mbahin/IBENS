#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import argparse
from scipy import stats
import numpy as np
import pandas as pd
import string
# Loading modules
import qPCR_functions
from qPCR_functions import cycles_col

###################################################################################################
# Functions

def max_second_derivative_smart(fx):
    ddf = [0, 0] + [fx[i + 2] - 2 * fx[i] + fx[i - 2] for i in range(2, len(fx) - 2)] + [0, 0]
    return np.argmax(ddf)


def modelled_fluo(row, equation):
    # Getting the fitting parameters except the first one, the min that we set to 0
    params = ["P_" + letter + "_final_" + equation.__name__ for letter in list(string.ascii_lowercase)[1:equation.__code__.co_argcount - 1]]
    args = [0] + row[params].values.tolist()
    # Running the fitting equation with these parameters
    return equation(np.arange(40), *args)


def logistic5p(x, minimum, maximum, IP, slope, f):
    try:
        return minimum + (maximum - minimum) / (np.power(1 + np.exp(slope * (x - IP)), f))
    except FloatingPointError:
        # If an error is raised, return a NaN vector so that it will not be used
        return np.repeat([np.nan], len(x))

### LRE functions ###
"""
def get_middle_pt(row):
    ''' Function to determine the first cycle after half of the maximum fluo is reached. '''
    Fmax = row[cycles_col_bump_corr].max()
    Fmin = row[cycles_col_bump_corr].min()
    for c in range(40):
        if row["F" + str(c).zfill(2) + "_bump_corr"] > Fmax - Fmin:
            return c
"""


def compute_Ec_val(row):
    ''' Compute efficiency values. '''
    res = []
    for i in range(1, 40):
        res.append(row["F" + str(i).zfill(2) + "_start0"] / row["F" + str(i - 1).zfill(2) + "_start0"] - 1)
    return res


def compute_deltaE_Emax_Fmax(row, cycle_start, cycle_end):
    ''' Compute the linear regression between fluo and efficiency values. '''
    slope, intercept, r_value, p_value, std_err = stats.linregress(list(row.loc[["F" + str(x).zfill(2) + "_start0" for x in range(cycle_start, cycle_end)]]),
                                                                   list(row.loc[["E" + str(x).zfill(2) for x in range(cycle_start, cycle_end)]]))
    return slope, intercept, intercept / (-1 * slope)


def compute_1_F0(row, c, Fmax, Emax):
    ''' Compute the F0 for one cycle. '''
    return Fmax / (1 + ((Fmax / row["F" + str(c).zfill(2) + "_start0"]) - 1) * np.power(Emax + 1, c))


def compute_F0(row):
    ''' Computes the final F0 for a curve. '''
    cycle_start = row["Rising_cycle"] + 1
    #cycle_start = row["Win_start"]
    for cycle_end in range(cycle_start + 4, 39):  ## Not 40 because the next cycle F0 will have to be tested
        # Compute deltaE (slope) and Emax (intercept) through linear regression (eq3) and Fmax (eq4)
        deltaE, Emax, Fmax = compute_deltaE_Emax_Fmax(row, cycle_start, cycle_end)
        F0s = []
        # Computing the F0 all cycles encompassed in the window + the following one
        for c in range(cycle_start, cycle_end + 1):
            F0s.append(compute_1_F0(row, c, Fmax, Emax))
        # Checking whether the last cycle is diverging from the mean of the others
        if np.abs(F0s[-1] / np.mean(F0s[:-1])) < 0.94:  # 0.94 = 1 - 0.06
            return deltaE, Emax, Fmax, cycle_end, cycle_end - cycle_start, np.median(F0s[:-1])
    return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
### End LRE functions ###

def get_new_rising_cycle(row):
    fluo = row[cycles_col_bump_corr]
    cycle = qPCR_functions.area_under_chord_pd(fluo)[1]
    f = ["F" + str(x).zfill(2) for x in range(cycle - 5, cycle + 5)]
    return max_second_derivative_smart(row[f]) + cycle - 5


def create_well_pairs(df):
    ''' Creates the appropriate well pairs to test the ratio compared to the ground truth. '''
    # Filtering out useless columns
    dataset = df.loc[:, ["Condition", "Gene", "Theoretical_N0_log", "Cq_method_N0_log", "Cy0", "F0_5p", "F0_LRE", "DL_predicted_N0_log", "key"]].copy()
    dataset.reset_index(inplace=True)
    # Auto-merging the dataset
    pairs = pd.merge(dataset, dataset, on="key", suffixes=["_1", "_2"])
    # Filtering out useless pairs  ## Creating the symmetrical pairs is useful when we want to make stats per genes for example
    pairs = pairs.loc[((pairs.Condition_1 != pairs.Condition_2) | (pairs.Gene_1 != pairs.Gene_2)) & (pairs.Theoretical_N0_log_1 != pairs.Theoretical_N0_log_2) & (pairs.WellID_1 > pairs.WellID_2)].copy()
    return pairs


def compute_t_test(df, gene1, condition1, gene2, condition2, feat, reversed=False):
    ''' Computes a t-test between the 6 replicates of two gene/conditoions to determine who is the greatest. '''
    # If we expect the greatest to have lesser values then we reverse the results.
    if reversed:
        rev = -1
    else:
        rev = 1
    # Collecting the two gene/condition values
    df1 = df.loc[(df.Gene == gene1) & (df.Condition == condition1), ["Gene", "Condition", "Theoretical_N0_log", feat]].copy()
    df2 = df.loc[(df.Gene == gene2) & (df.Condition == condition2), ["Gene", "Condition", "Theoretical_N0_log", feat]].copy()
    # Computing the T-test
    t, p = stats.ttest_ind(df1[feat].to_list(), df2[feat].to_list())
    # If the corrected p-value is greater than 0.05, the result is undetermined
    if ((p * gene_cond_pairs.shape[0]) / 2 > 0.05):
        return 0
    # Determining the result according to gene/condition ranking and T-test result (could probably be done in a smarter way...)
    elif df1["Theoretical_N0_log"].unique() > df2["Theoretical_N0_log"].unique():
        if t > 0:
            return 1 * rev
        else:
            return -1 * rev
    else:
        if t < 0:
            return 1 * rev
        else:
            return -1 * rev


def determine_log_cluster(row, lower_limit=0.2, upper_limit=2):
    ''' Function to determine the "log_cluster" (large, medium or small) of the comparison. '''
    log_diff = np.abs(row.Theoretical_N0_log_1 - row.Theoretical_N0_log_2)
    if log_diff >= upper_limit:
        return "large"
    elif log_diff >= lower_limit:
        return "intermediate"
    else:
        return "small"


def add_comparison(df, feat, method, ttest_res, reversed=False):
    ''' Function to add a method comparison. '''
    ## Comparison 1: computing the log ratio absolute diff
    # Adding the considered method ratio
    df.loc[:, method + "_ratio"] = np.log(df[feat + "_1"] / df[feat + "_2"])
    # Computing the difference with the theoretical ratio
    #df.loc[:, method + "_ratio_diff"] = df["Theoretical_N0_log_ratio"] - df[method + "_ratio"]
    # Computing the absolute difference with the theoretical ratio
    df.loc[:, method + "_abs_ratio_diff"] = np.abs(df["Theoretical_N0_log_ratio"] - df[method + "_ratio"])
    ## Comparison 2: expected pairwise ranking
    # Determining if the 2 wells ranking is true
    if reversed:
        df.loc[:, method + "_status"] = ((df.Theoretical_N0_log_1 > df.Theoretical_N0_log_2) & (df[feat + "_1"] < df[feat + "_2"])) | (
                    (df.Theoretical_N0_log_1 < df.Theoretical_N0_log_2) & (df[feat + "_1"] > df[feat + "_2"]))
    else:
        df.loc[:, method + "_status"] = ((df.Theoretical_N0_log_1 > df.Theoretical_N0_log_2) & (df[feat + "_1"] > df[feat + "_2"])) | (
                    (df.Theoretical_N0_log_1 < df.Theoretical_N0_log_2) & (df[feat + "_1"] < df[feat + "_2"]))
    ## Comparison 3: t-test between 6 replicates
    ttest_res_method = gene_cond_pairs.copy()
    ttest_res_method["Method"] = method
    ttest_res_method.loc[:, "6rep_test_res"] = gene_cond_pairs.apply(lambda row: compute_t_test(data, row.Gene_1, row.Condition_1, row.Gene_2, row.Condition_2, feat, reversed=reversed),                                                                   axis=1)
    ttest_res = pd.concat([ttest_res, ttest_res_method])
    # Displaying the results
    print("Method:", method)
    #print("Median ratio difference:", "{:.2f}".format(df[method + "_ratio_diff"].median()))
    print("Median ratio absolute difference:", "{:.2f}".format(df[method + "_abs_ratio_diff"].median()))
    print("Accuracy:", "{:.1%}".format(df.loc[df[method + "_status"], :].shape[0] / df.shape[0]))
    return df, ttest_res


def compute_Cy0(row):
    ''' Computes the Cy0 value. '''
    Richards_eq = "Richards"
    c = row["P_c_final_" + Richards_eq]
    b = row["P_d_final_" + Richards_eq]
    d = row["P_e_final_" + Richards_eq]
    Fb = row["P_a_final_" + Richards_eq]
    Fmax = row["P_b_final_" + Richards_eq]
    return c + b * np.log(d) - b * ((d + 1) / d) * (1 - ((Fb * np.power((d + 1) / d, d)) / Fmax))

###################################################################################################
# Main
if __name__ == "__main__":

    # Setting the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--plate", dest="plate", choices=["research", "QV"], help="The input plate to analyze.")
    options = parser.parse_args()

    # Setting data path according to the plate
    if options.plate == "research":
        input = "/data/biocomp/bahin/qPCR/Research_plate/Fixed_keys/final_data.total.REF.tsv"
        raw_data = pd.read_csv(input, sep="\t", index_col=0)
        data = raw_data.loc[(raw_data.Sigmoid_curve) & (raw_data.Pre_amplification == 0) & (raw_data.Gene != "Notch1a")
                            & ~((raw_data.Gene == "Crabp2a") & (raw_data.S_or_C == "PC"))
                            & ~((raw_data.Gene == "BMP4") & (raw_data.S_or_C == "PC")), :].copy()
        DL_input = "/data/biocomp/bahin/qPCR/Research_plate/Fixed_keys/DeepLearning/DL_predictions.research_plate.REF.tsv"
    else:  # QV plate
        input = "/data/biocomp/bahin/qPCR/Research_plate/Fixed_keys/QV_plate/final_data.total.REF.tsv"
        raw_data = pd.read_csv(input, sep="\t", index_col=0)
        data = raw_data.loc[raw_data.Sigmoid_curve & raw_data.Condition.isin(["PS" + str(x) for x in range(1, 13)] + ["PC" + str(x) for x in range(1, 13)]), :].copy()
        DL_input = "/data/biocomp/bahin/qPCR/Research_plate/Fixed_keys/DeepLearning/DL_predictions.QV_plate.REF.tsv"

    # Filtering out gene/condition without the 6 replicates showing a sigmoid curve
    data = data.groupby(["Gene", "Condition"]).filter(lambda row: row["Sample"].count() == 6)
    #data = data.sample(n=100, axis="index", random_state=42).copy()

    # Removing the bump ## TO MOVE TO qPCR_functions.py if validated with AG
    cycles_col_bump_corr = ["F" + str(i).zfill(2) + "_bump_corr" for i in range(40)]
    param_bump1 = "P_e_final_Gompertz_with_bump"
    param_bump2 = "P_f_final_Gompertz_with_bump"
    data[cycles_col_bump_corr] = qPCR_functions.correct_bump(data[cycles_col], data[param_bump1], data[param_bump2])

    # Computing the rising cycle new ## TO MOVE TO qPCR_functions.py if validated with AG
    data.loc[:, "Rising_cycle_new"] = data.apply(get_new_rising_cycle, axis=1)

    # Choosing the analyses to run
    HiFit = False
    #Cq_method, Cy0, fivep, LRE = False, False, False, False
    Cq_method, Cy0, fivep, LRE = True, True, True, True
    deepLearning = True

    # Creating cycle col names that start at 0 if logistic5p or LRE is processed
    if fivep or LRE:
        cycles_col_start0 = ["F" + str(i).zfill(2) + "_start0" for i in range(40)]

    ## Ct method
    if Cq_method:
        print("Running the Cq method...")
        """
        genes = ["Bactine2", "DeltaC", "DeltaD", "FGF8a", "Her7", "Hoxa1a", "Hoxa5a", "Hoxb1b", "Hoxd4a", "RPL13a1", "RPL13a2", "Rarab"]  # BMP4, Crapb2a, Notch1a and Cyp26a1 are discarded
        for ref_gene in genes:  # Developed to compare with different ref genes, could be simplified otherwise
            print("Ref_gene:", ref_gene)
            # Computing ref gene linear regression parameters
            slope, intercept, r_value, p_value, std_err = stats.linregress(data.loc[data.Gene == ref_gene, "Rising_cycle_new"], data.loc[data.Gene == ref_gene, "Theoretical_N0_log"])
            #  Computing predicted N0 (log)
            data.loc[data.Gene != ref_gene, "Cq_method_N0_log"] = data.loc[data.Gene != ref_gene, "Rising_cycle_new"] * slope + intercept
            # Creating the testing dataset (well pairs)  ## WARNING: comparing both ways?
            pairs = create_well_pairs(data.loc[data.Gene != ref_gene, :])
            #pairs.loc[:, ["Ref_gene", "Slope", "Intercept"]] = ref_gene, slope, intercept
            pairs = pairs.assign(Ref_gene=ref_gene, Slope=slope, Intercept=intercept)
            # Computing ratios
            pairs.loc[:, "Theoretical_N0_log_ratio"]  = np.log(pairs["Theoretical_N0_log_1"] / pairs["Theoretical_N0_log_2"])
            pairs.loc[:, "Cq_method_N0_log_ratio"] = np.log(pairs["Cq_method_N0_log_1"] / pairs["Cq_method_N0_log_2"])
            pairs.loc[:, "Cq_method_ratio_diff"] = np.abs(pairs.loc[:, "Theoretical_N0_log_ratio"] - pairs.loc[:, "Cq_method_N0_log_ratio"])
            # Accumulating results for each gene as the reference one
            if ref_gene == genes[0]:
                final = pairs.copy()
            else:
                final = pd.concat([final, pairs], join="inner")
        final.to_csv("/data/biocomp/bahin/qPCR/Research_plate/Fixed_keys/Cq_method.pairs.ratio_diff.tsv", sep="\t")
        """
        # Computing ref gene linear regression parameters
        ref_gene = "Bactine2"
        slope, intercept, r_value, p_value, std_err = stats.linregress(data.loc[data.Gene == ref_gene, "Rising_cycle_new"], data.loc[data.Gene == ref_gene, "Theoretical_N0_log"])
        #  Computing predicted N0 (log)
        data.loc[data.Gene != ref_gene, "Cq_method_N0_log"] = data.loc[data.Gene != ref_gene, "Rising_cycle_new"] * slope + intercept

    ## Cy0 method
    if Cy0:
        print("Running the Cy0 method...")
        # Computing Cy0 value
        data.loc[:, "Cy0"] = data.apply(compute_Cy0, axis=1)

    ## Logistic 5p method
    if fivep:
        print("Running logistic 5p method...")
        # Getting the modelled curve
        cycles_col_mod_5p = ["F" + str(i).zfill(2) + "_mod_5p" for i in range(40)]
        data[cycles_col_mod_5p] = data.apply(lambda row: modelled_fluo(row, logistic5p), axis=1, result_type="expand")
        # Computing the SDM
        data.loc[:, "5p_mod_SDM"] = data.apply(lambda row: max_second_derivative_smart(row[cycles_col_mod_5p]), axis=1)
        # Subtracting background
        data[cycles_col_start0] = data[cycles_col_bump_corr].sub(data["P_a_final_logistic5p"], axis=0)
        # Getting the efficiency at SDM
        data.loc[:, "Eff_5p_SDM"] = data.apply(lambda row: row["F" + str(row["5p_mod_SDM"]).zfill(2) + "_start0"] / row["F" + str(int(row["5p_mod_SDM"]) - 1).zfill(2) + "_start0"], axis=1)
        # Getting the F0 at SDM
        data.loc[:, "F0_5p"] = data.apply(lambda row: row["F" + str(row["5p_mod_SDM"]).zfill(2) + "_start0"] / np.power(row["Eff_5p_SDM"], row["5p_mod_SDM"]), axis=1)

    ## LRE method
    if LRE:
        print("Running LRE method...")
        # Baseline removal
        data[cycles_col_start0] = data[cycles_col_bump_corr].sub(data["P_a_final_Gompertz_with_bump"], axis=0)
        #data["Win_start"] = data.apply(get_middle_pt, axis=1)
        # Compute Ec values (eq5)
        Ec_cycles = ["E" + str(i).zfill(2) for i in range(1, 40)]
        data[Ec_cycles] = data.apply(compute_Ec_val, axis=1, result_type="expand")
        # Compute F0 through the recursive process
        data[["deltaE", "Emax_LRE", "Fmax_LRE", "LRE_win_end", "Win_size", "F0_LRE"]] = data.apply(compute_F0, axis=1, result_type="expand")
        #print(data.shape, data.dropna().shape)
        #print("Window mean size:", data["Win_size"].mean())

    ## HiFit method
    if HiFit:
        print("Running HiFit method...")

    ## Deep learning
    if deepLearning:
        print("Loading deep learning method data...")
        # Loading DL results
        dl = pd.read_csv(DL_input, sep="\t", index_col=0, header=None, names=["WellID", "DL_predicted_N0_log"])
        data = data.join(dl)

    # Creating all curves pairs
    data.to_csv("methods_results.tsv", sep="\t")
    print("Comparing methods...")
    #pairs = create_well_pairs(data)
    pairs = create_well_pairs(data.loc[data.Gene != "Bactine2", :])
    # Creating the theoretical N0 log ratio
    pairs.loc[:, "Theoretical_N0_log_ratio"] = np.log(pairs["Theoretical_N0_log_1"] / pairs["Theoretical_N0_log_2"])
    # Creating all valid gene/conditions pairs (either not the same gene or condition and only one comparison for 2 gene/condition, not both ways)
    ttest_res = pd.DataFrame()
    gene_cond_pairs = pairs.loc[((pairs.Gene_1 != pairs.Gene_2) | (pairs.Condition_1 != pairs.Condition_2)) & (pairs.Gene_1 > pairs.Gene_2), ["Gene_1", "Condition_1", "Theoretical_N0_log_1", "Gene_2", "Condition_2", "Theoretical_N0_log_2"]].drop_duplicates()
    gene_cond_pairs.loc[:, "log_cluster"] = gene_cond_pairs.apply(determine_log_cluster, axis=1)
    # Comparing Ct method (Warning, the well S55-A38, almost flat, gives a negative N0 log, raising an error in the function)
    if Cq_method:
        pairs, ttest_res = add_comparison(pairs.loc[(pairs.Gene_1 != ref_gene) & (pairs.Gene_2 != ref_gene), :].copy(), "Cq_method_N0_log", "Cq_method", ttest_res)
    # Comparing Cy0 method
    if Cy0:
        pairs, ttest_res = add_comparison(pairs, "Cy0", "Cy0_method", ttest_res, reversed=True)
    # Comparing logistic 5p method
    if fivep:
        pairs, ttest_res = add_comparison(pairs, "F0_5p", "5p_method", ttest_res)
    # Comparing LRE method
    if LRE:
        pairs, ttest_res = add_comparison(pairs, "F0_LRE", "LRE_method", ttest_res)
    # Comparing HiFit method
    if HiFit:
        pairs, ttest_res = add_comparison(pairs, "P_c_final_Gompertz", "HiFit_method", ttest_res, reversed=True)
    # Comparing the deep learning method
    if deepLearning:
        pairs, ttest_res = add_comparison(pairs, "DL_predicted_N0_log", "DeepLearning_method", ttest_res)
    pairs.to_csv("pairs.tsv", sep="\t")
    ttest_res["6rep_test_res"].replace({-1: "Wrong", 0: "Undetermined", 1: "Correct"}, inplace=True)
    ttest_res = ttest_res.groupby(["log_cluster", "Method"])["6rep_test_res"].value_counts(normalize=True).mul(100)
    ttest_res.to_csv("ttres.tsv", sep="\t", header=True)
