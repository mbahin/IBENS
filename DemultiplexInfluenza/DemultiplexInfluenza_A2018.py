#!/usr/bin/python

__author__ = "bahin"
""" Script to parse Read1 FASTQ from Influenza service. """

import argparse
from multiprocessing import Pool
from Bio import SeqIO
import numpy as np
import pandas as pd
import pysam
from natsort import natsorted
import matplotlib.pyplot as plt
import seaborn as sns
import sys


##### Functions #####
def load_FASTA_barcode(barcode_file):
    """ Load a FASTA file with one of the barcode sequences. """
    output = {}
    # Browsing the FASTA file
    for record in SeqIO.parse(barcode_file, "fasta"):
        output[str(record.seq)] = record.id
    return output

def get_RNA_fragment(sample_label, column):
    """ Get the mapped strain and RNA fragment info from a SAM file. """
    amplicons = {}
    # Loading the SAM file
    SAM_file = pysam.AlignmentFile(samples.at[sample_label, column], "rb")
    # Browsing the SAM file
    for read in SAM_file.fetch():
        # If read is unmapped
        if read.is_unmapped:
            amplicon = "Unmapped_Unmapped"
        # If read is not primary
        elif read.is_secondary:
            continue
        # If read is a multimapper
        elif read.query_name in amplicons:
            # If the previous mapping of this read was on a different strain and RNA fragment
            if amplicons[read.query_name] != read.reference_name:
                amplicon = "MappingAmbiguous_MappingAmbiguous"
            else:
                amplicon = read.reference_name
        else:
            amplicon = read.reference_name
        amplicons[read.query_name] = amplicon
    # Transforming the dict into a pandas DataFrame
    amplicons_df = pd.DataFrame.from_dict(amplicons, orient="index")  # We can rename the columns only from version 0.23.0
    amplicons_df.columns = ["Mapped"]
    # Splitting the info into strain and RNA fragment
    amplicons_df["Strain"], amplicons_df["RNA_fragment"] = amplicons_df["Mapped"].str.split("_").str
    amplicons_df.drop("Mapped", axis=1, inplace=True)
    return amplicons_df

def match_strain_info_from_R1_R2(row):
    """ Gather strain information retrieved from Read1 and Read2. """
    if (row.Strain_R1 == row.Strain_R2) and (row.Strain_R1 != "Unmapped"):
        return row.Strain_R1
    elif (row.Strain_R2 in ["", "Unmapped"]) and (row.Strain_R1 in ["H1N1", "H3N2"]):
        return row.Strain_R1
    elif (row.Strain_R1 in ["", "Unmapped"]) and (row.Strain_R2 in ["H1N1", "H3N2"]):
        return row.Strain_R2
    elif (row.Strain_R2 in ["", "MappingAmbiguous"]) and (row.Strain_R1 in ["H1N1", "H3N2"]):
        return row.Strain_R1
    elif (row.Strain_R1 in ["", "MappingAmbiguous"]) and (row.Strain_R2 in ["H1N1", "H3N2"]):
        return row.Strain_R2
    else:
        return "Undefined"

def match_RNA_fragment_info_from_R1_R2(row):
    """ Gather RNA fragment information retrieved from Read1 and Read2. """
    if (row.Strain_R1 == row.Strain_R2) and (row.Strain_R1 != "Unmapped"):
        if row.RNA_fragment_R1 == row.RNA_fragment_R2:
            return row.RNA_fragment_R1
        else:
            return "Undefined"
    elif (row.Strain_R2 in ["", "Unmapped"]) and (row.Strain_R1 in ["H1N1", "H3N2"]):
        return row.RNA_fragment_R1
    elif (row.Strain_R2 in ["", "Unmapped"]) and (row.Strain_R2 in ["H1N1", "H3N2"]):
        return row.RNA_fragment_R2
    elif (row.Strain_R2 in ["", "MappingAmbiguous"]) and (row.Strain_R1 in ["H1N1", "H3N2"]):
        return row.RNA_fragment_R1
    elif (row.Strain_R1 in ["", "MappingAmbiguous"]) and (row.Strain_R2 in ["H1N1", "H3N2"]):
        return row.RNA_fragment_R2
    elif (row.Strain_R1 != row.Strain_R2) and (row.RNA_fragment_R1 == row.RNA_fragment_R2):
        return row.RNA_fragment_R1
    else:
        return "Undefined"

# Browsing samples
def process_sample(sample_label):
    """ Full process one of the samples. """
    print("Processing " + sample_label + "...")
    print("Threshold for " + sample_label + ": " + str(samples.at[sample_label, "BC_threshold"]))
    # Defining the output directory
    output_dir = sample_label + "/"
    #output_dir = "Test/"

    ### Getting the barcodes from Read1
    barcodes = {}
    with open(output_dir + "/Indexes/" + sample_label + "_R1.3trimmed.76nt.fastq", "w") as output_file:
        # Browsing the 3' end trimmed Read1 FASTQ file
        print ("Barcode file: " + samples.at[sample_label, "Read1_Barcode"])
        for record in SeqIO.parse(samples.at[sample_label, "Read1_Barcode"], "fastq"):
            # Checking the trimmed read length (76nt is the exact expected length)
            if len(record.seq) == 76:
                SeqIO.write(record, output_file, "fastq")
                # Checking overhangB, overhangC and overhangD and registering the barcodes
                if (str(record.seq[16:20]) == "TTCG") and (str(record.seq[36:40] == "TGAC")) and (str(record.seq[56:60] == "ACCA")):
                    if (str(record.seq[0:16]) in bcA) and (str(record.seq[20:36]) in bcB) and (str(record.seq[40:56]) in bcC) and (str(record.seq[60:76]) in bcD):
                        barcodes[record.id[:-2]] = "A" + bcA[str(record.seq[0:16])] + "-B" + bcB[str(record.seq[20:36])] + "-C" + bcC[str(record.seq[40:56])] + "-D" + bcD[str(record.seq[60:76])]
    # Transforming the dict into a pandas DataFrame
    barcodes_df = pd.DataFrame.from_dict(barcodes, orient="index")
    barcodes_df.columns = ["Barcode"]

    ### Getting the UMI from Read1
    UMIs = {}
    # Browsing the 5' end trimmed Read1 FASTQ file
    print("UMI file: " + samples.at[sample_label, "Read1_UMI"])
    for record in SeqIO.parse(samples.at[sample_label, "Read1_UMI"], "fastq"):
        UMIs[record.id[:-2]] = str(record.seq)[0:9]
    # Transforming the dict into a pandas DataFrame
    UMIs_df = pd.DataFrame.from_dict(UMIs, orient="index")  # We can rename the columns only from version 0.23.0
    UMIs_df.columns = ["UMI"]

    ### Getting the amplicon from Read1 and Read2
    # Getting the info from Read1
    print("Amplicons file (from Read1): " + samples.at[sample_label, "Read1_RNA_fragment"])
    amplicons_R1_df = get_RNA_fragment(sample_label, "Read1_RNA_fragment")
    # Getting the info from Read2
    print("Amplicons file (from Read2): " + samples.at[sample_label, "Read2_RNA_fragment"])
    amplicons_R2_df = get_RNA_fragment(sample_label, "Read2_RNA_fragment")
    # Merging amplicon info obtained from Read1 and Read2
    amplicons_df = amplicons_R1_df.join(amplicons_R2_df, lsuffix="_R1", rsuffix="_R2", how="outer")
    # Replacing the NaN (from the outer join)
    amplicons_df.replace(np.nan, "", inplace=True)
    # Matching the strain and RNA fragment info from the 2 sources
    amplicons_df["Strain"] = amplicons_df.apply(match_strain_info_from_R1_R2, axis=1)
    amplicons_df["RNA_fragment"] = amplicons_df.apply(match_RNA_fragment_info_from_R1_R2, axis=1)
    # Reporting sample mapping statistics
    sample_res = amplicons_df.groupby(["Strain", "RNA_fragment"]).size()

    # Merging the 3 information
    reads = UMIs_df.join(amplicons_df).join(barcodes_df).sort_values(by=["Barcode", "UMI", "Strain", "RNA_fragment"])

    """
    # Filtering out barcode with too much UMIs
    if sample_label in ["i02", "i03"]:
        reads = reads.groupby("Barcode").filter(lambda x: x["UMI"].nunique() < 500)
    elif sample_label == "i04":
        reads = reads.groupby("Barcode").filter(lambda x: x["UMI"].nunique() < 1000)
    elif sample_label == "i05":
        reads = reads.groupby("Barcode").filter(lambda x: x["UMI"].nunique() < 1500)
    """

    ### Writing the output files
    #reads.set_index(["Index", "UMI", "Strain", "RNA_fragment"], inplace=True)  # Bad sorting, can't sort it like the others because the index "Barcode/UMI/Strain/RNA_fragment" is not unique
    #reads = reads.reindex(index=natsorted(reads.index))
    reads.to_csv(output_dir + sample_label + ".tsv", sep="\t")
    # Dropping old Read1 and Read2 columns
    reads.drop(["Strain_R1", "Strain_R2", "RNA_fragment_R1", "RNA_fragment_R2"], axis=1, inplace=True)
    # Computing the read counts before filtering
    droplet_read_counts = pd.DataFrame()
    droplet_read_counts["Total_reads_raw"] = reads.groupby("Barcode").UMI.count()
    droplet_read_counts["Total_unique_reads_raw"] = reads.groupby("Barcode").UMI.nunique()
    # Drawing the barplot to determine the appropriate threshold on the number of reads per barcode
    barplot = droplet_read_counts.reset_index()
    barplot = barplot.groupby("Total_unique_reads_raw", as_index=False).count()
    barplot["Weights"] = barplot["Barcode"] * barplot["Total_unique_reads_raw"]
    plt.hist(barplot["Total_unique_reads_raw"], np.logspace(0, 4, 100), weights=barplot["Weights"])
    plt.xscale("log")
    plt.xlabel("#Reads/barcode")
    plt.ylabel("count x #reads/BC")
    plt.title("Normalised reads count per barcode")
    plt.savefig(output_dir + sample_label + ".barplot.pdf", format="pdf")
    plt.clf()
    # Computing the droplet details (counts per barcode/UMI/strain/RNA_fragment) before filtering
    stats = pd.DataFrame()
    stats["Total_reads_count"] = reads.groupby(["Barcode", "UMI", "Strain", "RNA_fragment"]).UMI.count()
    """
    # Writing only the droplets with multiple strains
    final_noUnmapped = reads.loc[reads["Strain"] != "Unmapped"]
    multiStrains = pd.Series(final_noUnmapped.groupby(["Barcode"]).filter(lambda x: x.Strain.nunique() > 1)["Barcode"].unique())
    reads.loc[reads["Barcode"].isin(multiStrains)].groupby(["Barcode", "UMI", "Strain", "RNA_fragment"]).size().to_csv(options.sample_label + ".multi_strains.tsv", sep="\t")
    """
    # Filtering out fake barcodes (with very low read count)
    reads_BC_filter = reads.groupby(["Barcode"]).filter(lambda x: x["UMI"].count() >= samples.at[sample_label, "BC_threshold"])
    # Writing the filtered reads info
    reads_BC_filter.to_csv(output_dir + sample_label + ".all_filter.tsv", sep="\t")
    # Computing the reads counts after fake barcodes filtering
    droplet_read_counts["Total_reads_BC_filter"] = reads_BC_filter.groupby("Barcode").UMI.count()
    droplet_read_counts["Total_unique_reads_BC_filter"] = reads_BC_filter.groupby("Barcode").UMI.nunique()
    droplet_read_counts.reindex(index=natsorted(droplet_read_counts.index), copy=False)
    droplet_read_counts.to_csv(output_dir + sample_label + ".droplet_total_reads.tsv", sep="\t")
    # Computing the droplet details (counts per barcode/UMI/strain/RNA_fragment) after fake barcodes filtering
    stats["Total_reads_count_BC_filter"] = reads_BC_filter.groupby(["Barcode", "UMI", "Strain", "RNA_fragment"]).UMI.count()
    # Sorting the barcodes (natural sort)
    stats = stats.reindex(index=natsorted(stats.index))
    stats.to_csv(output_dir + sample_label + ".droplet_stats.tsv", sep="\t")
    # Drawing the clustermap to check reads per barcode filtering
    threshold_check_clustermap = stats.pivot_table(index="Barcode", columns=["Strain", "RNA_fragment"], values="Total_reads_count_BC_filter", aggfunc=np.sum, fill_value=0)
    ax = sns.clustermap(threshold_check_clustermap.corr().replace(np.nan, 0)).fig.suptitle("Amplicons within droplets correlation clustermap (sample " + sample_label + ")")
    plt.savefig(output_dir + sample_label + ".threshold_check_clustermap.pdf", format="pdf")
    plt.clf()
    # Drawing the clustermap of the discarded reads to check reads per barcode filtering
    reads_BC_discarded = reads.groupby(["Barcode"]).filter(lambda x: x["UMI"].count() < samples.at[sample_label, "BC_threshold"])
    stats_discarded = pd.DataFrame()
    stats_discarded["Total_reads_count_BC_filter"] = reads_BC_discarded.groupby(["Barcode", "UMI", "Strain", "RNA_fragment"]).UMI.count()
    threshold_check_clustermap_discarded = stats_discarded.pivot_table(index="Barcode", columns=["Strain", "RNA_fragment"], values="Total_reads_count_BC_filter", aggfunc=np.sum, fill_value=0)
    ax = sns.clustermap(threshold_check_clustermap_discarded.corr()).fig.suptitle("Amplicons within droplets correlation clustermap (sample " + sample_label + ")")
    plt.savefig(output_dir + sample_label + ".threshold_check_clustermap_discarded.pdf", format="pdf")
    plt.clf()
    # Drawing the scatter plot
    scatter_plot = reads_BC_filter.groupby(["Barcode", "Strain"], as_index=False).UMI.count()
    scatter_plot = scatter_plot.pivot(index="Barcode", columns="Strain", values="UMI")
    ax = sns.scatterplot(x="H1N1", y="H3N2", data=scatter_plot)
    plt.title("Strains associated to barcodes (sample " + sample_label + ")")
    plt.savefig(output_dir + sample_label + ".scatter_plot.pdf", format="pdf")
    plt.clf()
    # Counting the UMIs per barcode/strain/RNA_fragment
    UMI_counts = reads_BC_filter.groupby(["Barcode", "Strain", "RNA_fragment"]).count()
    # Adding a fake barcode (A97-B97-C97-D97) to be sure that the 8 fragments are present at least once for each strain
    UMI_counts.reset_index(inplace=True)
    for strain in ["H1N1", "H3N2"]:  # Improve the way this 16 fake rows are created
        for RNA_fragment in ["HA-1", "M-1", "NA-1", "NP-2", "NS-1", "PA-2", "PB1-1", "PB2-3"]:
            UMI_counts.loc[UMI_counts.shape[0]] = ["A97-B97-C97-D97", strain, RNA_fragment, 0]
    # Pivoting the DataFrame to produce the clustermap format
    clustermap = pd.pivot_table(UMI_counts, index="Barcode", columns=["Strain", "RNA_fragment"], fill_value=0)
    # Removing the fake barcode "A97-B97-C97-D97"
    clustermap.drop("A97-B97-C97-D97", inplace=True)
    # Writing the clustermap data
    clustermap.columns = clustermap.columns.droplevel()  # Dropping the UMI level that is common to every column
    clustermap.to_csv(output_dir + sample_label + ".clustermap.tsv", sep="\t")
    # Drawing and saving the clustermap
    sns.set(font_scale=0.7)
    ax = sns.clustermap(clustermap, col_cluster=False, vmax=25, cmap="magma_r", z_score=1)  # Drawing the clustermap
    plt.setp(ax.ax_heatmap.get_yticklabels(), rotation=0)  # Rotating the y labels
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=30)  # Rotating the x labels
    ax.ax_row_dendrogram.set_visible(False)  # Hiding the dendrogram
    plt.savefig(output_dir + sample_label + ".clustermap.pdf", format="pdf")
    plt.clf()
    # Creating the strains distribution within each RNA fragment
    strain_distribution = sample_res.reset_index()  # Resetting the index to be able to filter on the strain and RNA_fragment
    strain_distribution = strain_distribution.loc[(strain_distribution["Strain"] != "Undefined") & (strain_distribution["Strain"] != "MappingAmbiguous") & (strain_distribution["RNA_fragment"] != "Undefined") & (strain_distribution["RNA_fragment"] != "MappingAmbiguous")]
    strain_distribution.set_index(["Strain", "RNA_fragment"], inplace=True)  # Re-setting the index to be able to do a groupby
    strain_distribution = strain_distribution.groupby("RNA_fragment").apply(lambda x: 100 * x / float(x.sum()))
    strain_distribution.columns = ["Percentages"]
    strain_distribution.reset_index(inplace=True)  # Resetting the index to be able to pivot
    strain_distribution = strain_distribution.pivot(index="RNA_fragment", columns="Strain", values="Percentages")
    strain_distribution.plot.bar(stacked=True)
    plt.title("Strains_distribution_per_RNA_fragment (sample " + sample_label + ")")
    plt.savefig(output_dir + sample_label + ".strains_distribution_per_RNA_fragment.pdf", format="pdf")
    # Reporting the end of the analysis for this sample
    print(sample_label + " processed.")
    return sample_label, sample_res

##### End functions #####

if __name__ == "__main__":
    # Getting command-line options back
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", dest="samples", required=True, help="Samples config file")
    parser.add_argument("-a", dest="bcA", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/A_forward.fasta", help="BarcodeA FASTA file")
    parser.add_argument("-b", dest="bcB", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/B_forward.fasta", help="BarcodeB FASTA file")
    parser.add_argument("-c", dest="bcC", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/C_forward.fasta", help="BarcodeC FASTA file")
    parser.add_argument("-d", dest="bcD", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/D_forward.fasta", help="BarcodeD FASTA file")
    parser.add_argument("-t", dest="nb_threads", type=int, default=1, help="Number of threads to process the data")
    options = parser.parse_args()

    # Parsing the samples config file
    samples = pd.read_csv(options.samples, sep="\t", header=0, index_col=0)

    # Loading the barcodes
    bcA = load_FASTA_barcode(options.bcA)
    bcB = load_FASTA_barcode(options.bcB)
    bcC = load_FASTA_barcode(options.bcC)
    bcD = load_FASTA_barcode(options.bcD)

    # Run the analysis for each sample in parallel if more than one thread was called
    if options.nb_threads > 1:
        pool = Pool(options.nb_threads)
        results = list(pool.imap_unordered(process_sample, samples.index))
    else:  # Non-multithread version
        results = []
        for sample in samples.index:
            results.append(process_sample(sample))

    # Writing the file with all mapping statistics
    all_samples_mapping = pd.concat([r[1] for r in results], axis=1)
    all_samples_mapping.columns = [r[0] for r in results]
    all_samples_mapping.sort_index(axis=1, inplace=True)
    all_samples_mapping.to_csv("all_mapping_stats.tsv", sep="\t")

    # Producing the summary strains barplot
    #all_samples_mapping.groupby(["Strain"]).sum().T.plot(kind="bar", stacked=True, colormap="Set2", legend="reverse").legend(bbox_to_anchor=(1.5, 1))  # Not normalized
    all_samples_mapping.groupby(["Strain"]).sum().apply(lambda x: 100 * x / x.sum()).T.plot(kind="bar", stacked=True, colormap="Set2", legend="reverse").legend(bbox_to_anchor=(1.5, 1))
    plt.title("Strains distribution within samples")
    plt.savefig("strains_barplot.pdf", bbox_inches="tight")
    plt.clf()
    # Producing the amplicons barplot
    all_samples_mapping.reset_index(inplace=True)
    amplicons_barplot = all_samples_mapping.loc[(all_samples_mapping["Strain"] != "Undefined") & (all_samples_mapping["Strain"] != "MappingAmbiguous") & (all_samples_mapping["RNA_fragment"] != "Undefined") & (all_samples_mapping["RNA_fragment"] != "MappingAmbiguous")]
    amplicons_barplot = pd.melt(amplicons_barplot, id_vars=["Strain", "RNA_fragment"], value_vars=["i02", "i03", "i04", "i05", "i06", "i07", "i08"], var_name="Sample")
    amplicons_barplot["Amplicon"] = amplicons_barplot["Strain"] + "_" + amplicons_barplot["RNA_fragment"]
    ax = sns.barplot(data=amplicons_barplot, x="Amplicon", y="value", hue="Sample", ci=None)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title("Amplicons distribution within samples")
    plt.savefig("amplicons_barplot.pdf")
    plt.clf()
