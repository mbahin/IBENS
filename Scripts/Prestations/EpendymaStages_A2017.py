#!/usr/bin/python

__author__ = "bahin"
""" Scripts to produce files for EpendymaStages prestation. """

import argparse
import pandas as pd
import numpy as np
import collections

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument("-e", dest="EpenDiff", default="/projects/biocompan/Finished/201708_IBENS_Spassky_EpendymaStages_A2017/Bioinformatics_platform/Ependiff_A2016/Eoulsan/diffanaresultsannotation_deseq2_Experiment1-diffana_DIV4_vs_IP-PS6.reverse.tsv", help="EpenDiff Eoulsan output file")
parser.add_argument("-f", dest="FateChoice", default="/projects/biocompan/Finished/201708_IBENS_Spassky_EpendymaStages_A2017/Bioinformatics_platform/FateChoice_C2013/Eoulsan/diffanaresultsannotation_deseq2_Experiment1-diffana_DIV4_vs_DIV-2.tsv", help="FateChoice Eoulsan output file")
parser.add_argument("-g", dest="GemC1", default="/projects/biocompan/Bioinformatics_services/201904_IBENS_Spassky_EpendymaStages_A2019/Bioinformatics_platform/eoulsan/diffanaresultsannotation_output/diffanaresultsannotation_deseq2_Experiment1-diffana_GemC1_vs_H2B-GFP.tsv", help="GemC1 Eoulsan output file")
parser.add_argument("-m", dest="MCIDAS", default="/projects/biocompan/Bioinformatics_services/201904_IBENS_Spassky_EpendymaStages_A2019/Bioinformatics_platform/eoulsan/diffanaresultsannotation_output/diffanaresultsannotation_deseq2_Experiment1-diffana_MCIDAS_vs_H2B-GFP.tsv", help="MCIDAS Eoulsan output file")
parser.add_argument("-o", dest="output", required=True, help="Output filename")
options = parser.parse_args()

# Defininf genes lists
genes_lists_path = "/projects/biocompan/Finished/201708_IBENS_Spassky_EpendymaStages_A2017/Bioinformatics_platform/Genes_list/"
genes_lists = collections.OrderedDict([("G1_S_Dominguez2016", "G1_S_Dominguez2016.gene_names.formatted.tsv"),
                                       ("G2_M_Dominguez2016", "G2_M_Dominguez2016.gene_names.formatted.tsv"),
                                       ("Adult_neural_stem_cell_TableS1Codega2014.activated_NSC", "Adult_neural_stem_cell_TableS1Codega2014.activated_NSC.formatted.tsv"),
                                       ("Adult_neural_stem_cell_TableS1Codega2014.quiescent_NSC", "Adult_neural_stem_cell_TableS1Codega2014.quiescent_NSC.formatted.tsv"),
                                       ("TOPmRNA_Thoreen2012", "TOPmRNA_Thoreen2012.formatted.tsv"),
                                       ("Ftdye_Ng-type_Telley", "Nikita/Ftdye_Ng-type_Telley.formatted.tsv"),
                                       ("Ftdye_N-type_Telley", "Nikita/Ftdye_N-type_Telley.formatted.tsv"),
                                       ("Ftdye_P-type_Telley", "Nikita/Ftdye_P-type_Telley.formatted.tsv"),
                                       ("GPCR_Codega", "Nikita/GPCR_Codega.tsv"),
                                       ("Marker_Codega", "Nikita/Marker_Codega.formatted.tsv"),
                                       ("Prolif_Codega", "Nikita/Prolif_Codega.tsv"),
                                       ("PS6RP_Knight_2012", "Nikita/PS6RP_Knight_2012.formatted.tsv"),
                                       ("Quiescence_Codega", "Nikita/Quiescence_Codega.tsv"),
                                       ("TF_Codega", "Nikita/TF_Codega.formatted.tsv"),
                                       ("astrocyte.Zeisel", "Zeisel/astrocyte.Zeisel.txt"),
                                       ("CA1_Pyramidal.Zeisel", "Zeisel/CA1_Pyramidal.Zeisel.txt"),
                                       ("endothelial.Zeisel", "Zeisel/endothelial.Zeisel.txt"),
                                       ("ependymal.Zeisel", "Zeisel/ependymal.Zeisel.txt"),
                                       ("interneuron.Zeisel", "Zeisel/interneuron.Zeisel.txt"),
                                       ("microglia.Zeisel", "Zeisel/microglia.Zeisel.txt"),
                                       ("mural.Zeisel.txt", "Zeisel/mural.Zeisel.txt"),
                                       ("oligodendrocyte.Zeisel", "Zeisel/oligodendrocyte.Zeisel.txt"),
                                       ("S1_Pyramidal.Zeisel", "Zeisel/S1_Pyramidal.Zeisel.txt"),
                                       ("Senescence.Crowe", "senescence.Mus.gene_names.uniq.Crowe_2016.tsv"),
                                       ("Senescence.Hernandez", "senescence.Mus.gene_name.uniq.Hernandez-Segura-2017.tsv"),
                                       ("Senescence.Hernandez-Segura_RS", "senescence.Mus.gene_name.uniq.Hernandez-Segura_RS.tsv"),
                                       ("Senescence.Hernandez-Segura_OIS", "senescence.Mus.gene_name.uniq.Hernandez-Segura_OIS.tsv"),
                                       ("Senescence.Hernandez-Segura_IRIS", "senescence.Mus.gene_name.uniq.Hernandez-Segura_IRIS.tsv")])

# Parsing EpenDiff input file
data = pd.read_csv(options.EpenDiff, index_col=0, usecols=[0, 2, 6, 9], header=0, names=["Gene_ID", "EpenDiff_log2FC", "EpenDiff_adj_pvalue", "Gene_name"], sep="\t")
data = data[["Gene_name", "EpenDiff_log2FC", "EpenDiff_adj_pvalue"]]
data.loc[:, "EpenDiff_FC"] = 2 ** data["EpenDiff_log2FC"]
data.loc[data.EpenDiff_log2FC == 0, ["EpenDiff_FC", "EpenDiff_log2FC"]] = [np.nan] * 2  # Setting the "regular" and log2 FC to NaN if the log2FC was 0

# Parsing FateChoice input file
FateChoice = pd.read_csv(options.FateChoice, index_col=0, usecols=[0, 2, 6], header=0, names=["Gene_ID", "FateChoice_log2FC", "FateChoice_adj_pvalue"], sep="\t")
FateChoice.loc[:, "FateChoice_FC"] = 2 ** FateChoice["FateChoice_log2FC"]

# Parsing GemC1 input file
GemC1 = pd.read_csv(options.GemC1, index_col=0, usecols=[0, 2, 6], header=0, names=["Gene_ID", "GemC1_log2FC", "GemC1_adj_pvalue"], sep="\t")
GemC1.loc[:, "GemC1_FC"] = 2 ** GemC1["GemC1_log2FC"]

# Parsing MCIDAS input file
MCIDAS = pd.read_csv(options.MCIDAS, index_col=0, usecols=[0, 2, 6], header=0, names=["Gene_ID", "MCIDAS_log2FC", "MCIDAS_adj_pvalue"], sep="\t")
MCIDAS.loc[:, "MCIDAS_FC"] = 2 ** MCIDAS["MCIDAS_log2FC"]

# Merging input files
data = data.join([FateChoice, GemC1, MCIDAS])
data.drop(["EpenDiff_log2FC", "FateChoice_log2FC", "GemC1_log2FC", "MCIDAS_log2FC"], axis=1, inplace=True)

# Filtering only significant genes (for FateChoice and EpenDiff only)
EpenDiff_sig = data.loc[(data.EpenDiff_adj_pvalue < 0.01) & (data.EpenDiff_FC > 1.5), "Gene_name"].drop_duplicates()
data.loc[0:EpenDiff_sig.shape[0], "EpenDiff_significant_genes"] = EpenDiff_sig.values
FateChoice_sig = data.loc[(data.FateChoice_adj_pvalue < 0.01) & (data.FateChoice_FC > 2), "Gene_name"].drop_duplicates()
data.loc[0:FateChoice_sig.shape[0], "FateChoice_significant_genes"] = FateChoice_sig.values

# Browsing gene sets to order the genes correctly (only significant in FateChoice, only in EpenDiff, both and none)
for list_name in genes_lists:
    print(list_name)
    gene_set = pd.read_csv(genes_lists_path + genes_lists[list_name], header=None, names=[list_name])
    s = gene_set.merge(data, left_on=list_name, right_on="Gene_name", how="left")  # Left merge to get all the genes from the gene set
    s.loc[:, "Gene_name"] = s.loc[:, list_name]  # Copying gene names from the gene set that are not in the main DataFrame
    s1 = s.loc[(s.Gene_name.str.lower().isin(FateChoice_sig.str.lower())) & (~s.Gene_name.str.lower().isin(EpenDiff_sig.str.lower())), "Gene_name"].drop_duplicates()
    s2 = s.loc[(~s.Gene_name.str.lower().isin(FateChoice_sig.str.lower())) & (s.Gene_name.str.lower().isin(EpenDiff_sig.str.lower())), "Gene_name"].drop_duplicates()
    s3 = s.loc[(s.Gene_name.str.lower().isin(FateChoice_sig.str.lower())) & (s.Gene_name.str.lower().isin(EpenDiff_sig.str.lower())), "Gene_name"].drop_duplicates()
    s4 = s.loc[(~s.Gene_name.str.lower().isin(FateChoice_sig.str.lower())) & (~s.Gene_name.str.lower().isin(EpenDiff_sig.str.lower())), "Gene_name"].drop_duplicates()
    s = pd.concat([s1, s2, s3, s4])
    s.columns = list_name
    data.loc[0:s.shape[0], list_name] = s.values

# Writing the output file
data.to_csv(options.output, sep="\t")
"""
for list_name in genes_lists:
    # For the True/False columns (volcano plot)
    s = df_bool.loc[:, "Name"].isin(pd.Series(read_list(genes_lists_path + genes_lists[list_name])))
    s.name = list_name
    df_bool = pd.concat([df_bool, s], axis=1)
df_bool.to_csv(options.output + "2", sep="\t", index=False)
"""