#!/usr/bin/python

__author__ = "bahin"
""" Scripts to do a randomization test. """

import pandas as pd
import os
from scipy.stats import norm
import math

# Loading the fold changes
#DE_res_filepath = "/projects/biocompan/Finished/201708_IBENS_Spassky_EpendymaStages_A2017/Bioinformatics_platform/FateChoice_C2013/Eoulsan/diffanaresultsannotation_deseq2_Experiment1-diffana_DIV4_vs_DIV-2.tsv"
DE_res_filepath = "/projects/biocompan/Finished/201708_IBENS_Spassky_EpendymaStages_A2017/Bioinformatics_platform/Ependiff_A2016/Eoulsan/diffanaresultsannotation_deseq2_Experiment1-diffana_DIV4_vs_IP-PS6.reverse.tsv"
DE_res = pd.read_csv(DE_res_filepath, sep="\t", usecols=[0, 1, 2, 9])
DE_res.columns = ["Id", "baseMean", "FC", "Gene_name"]
# Lowering gene names to be comparable
DE_res["Gene_name"] = DE_res["Gene_name"].str.lower()
# Filtering out lowly expressed genes
DE_res = DE_res.loc[DE_res.baseMean > 5]

# Loading the gene sets
dirpath = "/projects/biocompan/Finished/201708_IBENS_Spassky_EpendymaStages_A2017/Bioinformatics_platform/Genes_list/GMT_format/"
gene_sets = os.listdir(dirpath)

# Getting gene set median fold change for genes in the DE results
res = pd.DataFrame()
for gene_set_filepath in gene_sets:
    #print(gene_set_filepath)
    # Loading gene set
    gene_set = pd.read_csv(dirpath + gene_set_filepath, header=None, sep="\t").T
    gene_set = gene_set.loc[2:, ].dropna()
    gene_set.columns = ["Gene_name"]
    # Lowering gene names to be comparable
    gene_set["Gene_name"] = gene_set["Gene_name"].str.lower()
    # Getting genes from the gene set present in the RNA-Seq file
    gene_set = gene_set.merge(DE_res, on="Gene_name")
    # Computing metrics for the CDF and SF
    gene_set_mean = gene_set["FC"].mean()
    DE_res_mean = DE_res["FC"].mean()
    DE_res_sigma = DE_res["FC"].std()
    gene_set_size = gene_set.shape[0]
    # Computing CDF (for up-regulated genes) and SF (for down-regulated genes)
    res.loc[gene_set_filepath, "CDF"] = norm().cdf((gene_set_mean - DE_res_mean) / (DE_res_sigma / math.sqrt(gene_set_size)))
    res.loc[gene_set_filepath, "SF"] = norm().sf((gene_set_mean - DE_res_mean) / (DE_res_sigma / math.sqrt(gene_set_size)))

# Correcting the p-values
res.loc[:, "CDF_corrected"] = res["CDF"] * len(gene_sets)
res.loc[:, "SF_corrected"] = res["SF"] * len(gene_sets)

# Storing the results
#res.to_csv("", sep="\t")