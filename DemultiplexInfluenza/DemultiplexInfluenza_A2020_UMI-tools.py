#!/usr/bin/python

__author__ = "bahin"
''' Script to parse data from Influenza service. '''

import argparse
from multiprocessing import Pool
from Bio import SeqIO
import numpy as np
import pandas as pd
import pysam
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import sys


##### Functions #####
def load_FASTA_barcode(barcode_file):
    ''' Load a FASTA file with one of the barcode sequences. '''
    output = {}
    # Browsing the FASTA file
    for record in SeqIO.parse(barcode_file, "fasta"):
        output[str(record.seq)] = record.id
    return output

def get_RNA_fragment(sample_label, column):
    ''' Get the mapped strain and RNA fragment info from a SAM file. '''
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
    ''' Gather strain information retrieved from Read1 and Read2. '''
    if (row.Strain_R1 == row.Strain_R2) and (row.Strain_R1 not in ["", "Unmapped", "MappingAmbiguous"]):
        return row.Strain_R1
    #elif (row.Strain_R2 in ["", "Unmapped"]) and (row.Strain_R1 in ["H1N1", "H3N2"]):
    elif (row.Strain_R2 in ["", "Unmapped", "MappingAmbiguous"]) and (row.Strain_R1 not in ["", "Unmapped", "MappingAmbiguous"]):
        return row.Strain_R1
    #elif (row.Strain_R1 in ["", "Unmapped"]) and (row.Strain_R2 in ["H1N1", "H3N2"]):
    elif (row.Strain_R1 in ["", "Unmapped", "MappingAmbiguous"]) and (row.Strain_R2 not in ["", "Unmapped", "MappingAmbiguous"]):
        return row.Strain_R2
    else:
        return "Undefined"

def match_RNA_fragment_info_from_R1_R2(row):
    ''' Gather RNA fragment information retrieved from Read1 and Read2. '''
    if (row.RNA_fragment_R1 == row.RNA_fragment_R2) and (row.RNA_fragment_R1 not in ["", "Unmapped", "MappingAmbiguous"]):
        return row.RNA_fragment_R1
    elif (row.RNA_fragment_R2 in ["", "Unmapped", "MappingAmbiguous"]) and (row.RNA_fragment_R1 not in ["", "Unmapped", "MappingAmbiguous"]):
        return row.RNA_fragment_R1
    elif (row.RNA_fragment_R1 in ["", "Unmapped", "MappingAmbiguous"]) and (row.RNA_fragment_R2 not in ["", "Unmapped", "MappingAmbiguous"]):
        return row.RNA_fragment_R2
    else:
        return "Undefined"

def weighted_hist(x, weights, **kwargs):
    ''' Wrapper for seaborn distplot function in order to pass the weights. '''
    sns.distplot(x, hist_kws={"weights": weights}, **kwargs)

# Browsing samples
def process_sample(sample_label):
    ''' Full process one of the samples. '''
    print("Processing " + sample_label + "...")
    print("Threshold for " + sample_label + ": " + str(samples.at[sample_label, "BC_threshold"]))
    # Defining the output directory
    output_dir = sample_label + "/"
    #output_dir = "Test/"

    ### Getting the barcodes from Read1
    barcodes = {}
    """
    with open(output_dir + "/Indexes/" + sample_label + "_R1.3trimmed.84nt.fastq", "w") as output_file:
        # Browsing the 3' end trimmed Read1 FASTQ file
        print ("Barcode file: " + samples.at[sample_label, "Read1_Barcode"])
        for record in SeqIO.parse(samples.at[sample_label, "Read1_Barcode"], "fastq"):
            # Checking the trimmed read length (76nt is the exact expected length)
            if len(record.seq) == 84:
                SeqIO.write(record, output_file, "fastq")
                # Checking overhangB, overhangC and overhangD and registering the barcodes
                if (str(record.seq[20:24]) == "TGAC") and (str(record.seq[40:44] == "ACCA")) and (str(record.seq[60:64] == "CAAC")):
                    if (str(record.seq[4:20]) in bcB) and (str(record.seq[24:40]) in bcC) and (str(record.seq[44:60]) in bcD) and (str(record.seq[64:84]) in bcE):
                        barcodes[record.id] = "B" + bcB[str(record.seq[4:20])] + "-C" + bcC[str(record.seq[24:40])] + "-D" + bcD[str(record.seq[44:60])] + "-E" + bcE[str(record.seq[64:84])]
    """

    def gather_BC(row):
        if np.isnan(row.B):
            B = "B??"
        else:
            B = "B" + str(int(row.B))
        if np.isnan(row.C):
            C = "C??"
        else:
            C = "C" + str(int(row.C))
        if np.isnan(row.D):
            D = "D??"
        else:
            D = "D" + str(int(row.D))
        if np.isnan(row.E):
            E = "E??"
        else:
            E = "E" + str(int(row.E))
        BC = B + "-" + C + "-" + D + "-" + E
        if np.isnan(row.B) or np.isnan(row.C) or np.isnan(row.D) or np.isnan(row.E):
            return "Undefined (" + BC + ")"
        else:
            return B + "-" + C + "-" + D + "-" + E


    filenameB = "/projects/biocompan/Bioinformatics_services/202003_Pasteur_Isel_DemultiplexInfluenza_A2020/Runs/MiSeq/UMI-tools/i06.BBC.fastq"
    filenameC = "/projects/biocompan/Bioinformatics_services/202003_Pasteur_Isel_DemultiplexInfluenza_A2020/Runs/MiSeq/UMI-tools/i06.CBC.fastq"
    filenameD = "/projects/biocompan/Bioinformatics_services/202003_Pasteur_Isel_DemultiplexInfluenza_A2020/Runs/MiSeq/UMI-tools/i06.DBC.fastq"
    filenameE = "/projects/biocompan/Bioinformatics_services/202003_Pasteur_Isel_DemultiplexInfluenza_A2020/Runs/MiSeq/UMI-tools/i06.EBC.fastq"
    for record in SeqIO.parse(filenameB, "fastq"):
        barcodes[record.id.split("_")[0]] = {}
        barcodes[record.id.split("_")[0]]["B"] = bcB[str(record.id.split("_")[1])]
    for record in SeqIO.parse(filenameC, "fastq"):
        if record.id.split("_")[0] not in barcodes:
            barcodes[record.id.split("_")[0]] = {}
        barcodes[record.id.split("_")[0]]["C"] = bcC[str(record.id.split("_")[1])]
    for record in SeqIO.parse(filenameD, "fastq"):
        if record.id.split("_")[0] not in barcodes:
            barcodes[record.id.split("_")[0]] = {}
        barcodes[record.id.split("_")[0]]["D"] = bcD[str(record.id.split("_")[1])]
    for record in SeqIO.parse(filenameE, "fastq"):
        if record.id.split("_")[0] not in barcodes:
            barcodes[record.id.split("_")[0]] = {}
        barcodes[record.id.split("_")[0]]["E"] = bcE[str(record.id.split("_")[1])]
    # Transforming the dict into a pandas DataFrame
    barcodes_df = pd.DataFrame.from_dict(barcodes, orient="index")
    barcodes_df = barcodes_df.apply(pd.to_numeric)
    barcodes_df["Barcode"] = barcodes_df.apply(gather_BC, axis=1)
    barcodes_df.to_csv(sample_label + ".BC.tsv", sep="\t")
    barcodes_df.drop(columns=["B", "C", "D", "E"], inplace=True)
    #barcodes_df.columns = ["Barcode"]

    ### Getting the UMI from Read1
    UMIs = {}
    # Browsing the 5' end trimmed Read1 FASTQ file
    print("UMI file: " + samples.at[sample_label, "Read1_UMI"])
    for record in SeqIO.parse(samples.at[sample_label, "Read1_UMI"], "fastq"):
        UMIs[record.id] = str(record.seq)[0:9]
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
    #amplicons_df.replace(np.nan, "", inplace=True)
    amplicons_df.replace(np.nan, "Undefined", inplace=True)
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
    print("Sample " + sample_label + " - Considered reads count: " + str(reads.shape[0]))
    print("Sample " + sample_label + " - Considered BC count: " + str(len(reads.groupby("Barcode"))))
    reads_BC_filter = reads.groupby(["Barcode"]).filter(lambda x: x["UMI"].count() >= samples.at[sample_label, "BC_threshold"])
    print("Sample " + sample_label + " - Considered reads count after barcode filtering: " + str(reads_BC_filter.shape[0]))
    print("Sample " + sample_label + " - Considered BC count after filtering: " + str(len(reads_BC_filter.groupby("Barcode"))))
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
    if not reads_BC_discarded.empty:  # Happens mostly if the BC threshold is set to 1
        stats_discarded = pd.DataFrame()
        stats_discarded["Total_reads_count_BC_filter"] = reads_BC_discarded.groupby(["Barcode", "UMI", "Strain", "RNA_fragment"]).UMI.count()
        threshold_check_clustermap_discarded = stats_discarded.pivot_table(index="Barcode", columns=["Strain", "RNA_fragment"], values="Total_reads_count_BC_filter", aggfunc=np.sum, fill_value=0)
        ax = sns.clustermap(threshold_check_clustermap_discarded.corr()).fig.suptitle("Amplicons within droplets correlation clustermap (sample " + sample_label + ")")
        plt.savefig(output_dir + sample_label + ".threshold_check_clustermap_discarded.pdf", format="pdf")
        plt.clf()
    else:
        print("Oh no, empty!")
    # Counting the UMIs per barcode/strain/RNA_fragment
    UMI_counts = reads_BC_filter.groupby(["Barcode", "Strain", "RNA_fragment"]).count()
    # Adding a fake barcode (A97-B97-C97-D97) to be sure that the 8 fragments are present at least once for each strain
    UMI_counts.reset_index(inplace=True)
    for strain in ["H1N1", "H3N2"]:  # Improve the way this 16 fake rows are created
        for RNA_fragment in RNA_fragments:
            UMI_counts.loc[UMI_counts.shape[0]] = ["B97-C97-D97-E97", strain, RNA_fragment, 0]
    # Filtering out undefined data (strain and/or RNA fragment)
    UMI_counts = UMI_counts.loc[(UMI_counts.Strain != "Undefined") & (UMI_counts.RNA_fragment != "Undefined"), :]
    # Pivoting the DataFrame to produce the clustermap format
    clustermap = pd.pivot_table(UMI_counts, index="Barcode", columns=["Strain", "RNA_fragment"], fill_value=0)
    # Removing the fake barcode "A97-B97-C97-D97"
    clustermap.drop("B97-C97-D97-E97", inplace=True)
    # Writing the clustermap data
    clustermap.columns = clustermap.columns.droplevel()  # Dropping the UMI level that is common to every column
    clustermap.to_csv(output_dir + sample_label + ".clustermap.tsv", sep="\t")
    clustermap = clustermap.loc[:, clustermap.sum(axis=0) != 0]  # Removing the columns that sum to 0 (FloatingPointError otherwise)
    # Normalizing the data by row and by col
    clustermap = clustermap.apply(lambda x: x/x.max(), axis=0)
    clustermap = clustermap.apply(lambda x: x/x.max(), axis=1)
    # Drawing and saving the clustermap
    sys.setrecursionlimit(10000)
    sns.set(font_scale=0.7)
    """
    ax = sns.clustermap(clustermap, col_cluster=False, cmap="magma_r")  # Drawing the clustermap
    plt.setp(ax.ax_heatmap.get_yticklabels(), rotation=0)  # Rotating the y labels
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=30)  # Rotating the x labels
    ax.ax_row_dendrogram.set_visible(False)  # Hiding the dendrogram
    plt.savefig(output_dir + sample_label + ".clustermap.pdf", format="pdf")
    plt.clf()
    sys.setrecursionlimit(1000)
    """

    # Reporting the end of the analysis for this sample
    print(sample_label + " processed.")
    return sample_label, sample_res

##### End functions #####

if __name__ == "__main__":
    # Getting command-line options back
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", dest="samples", required=True, help="Samples config file")
    parser.add_argument("-b", dest="bcB", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/B_forward.fasta", help="BarcodeB FASTA file")
    parser.add_argument("-c", dest="bcC", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/C_forward.fasta", help="BarcodeC FASTA file")
    parser.add_argument("-d", dest="bcD", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/D_forward.fasta", help="BarcodeD FASTA file")
    parser.add_argument("-e", dest="bcE", default="/projects/biocompan/Bioinformatics_services/202003_Pasteur_Isel_DemultiplexInfluenza_A2020/Data/E_forward.fasta", help="BarcodeD FASTA file")
    parser.add_argument("-t", dest="nb_threads", type=int, default=1, help="Number of threads to process the data")
    options = parser.parse_args()

    # Setting variables
    RNA_fragments = ["HA", "M", "NAseg", "NP", "NS", "PA", "PB1", "PB2"]
    seg_colors = {"H1N1": "#109CF8", "H3N2": "#EC8D79", "Undefined": "#B8B8BB"}
    colors = pd.Series([seg_colors["H1N1"], seg_colors["H3N2"], seg_colors["Undefined"]], index=["H1N1", "H3N2", "Undefined"])
    leg = [mpatches.Patch(color=seg_colors[c], label=c) for c in seg_colors]
    sample_loc = {"i04": [0, 0], "i05": [0, 1], "i06": [1, 0], "i08": [1, 1]}  ## Placing samples on subplots, not ideal at all...

    # Parsing the samples config file
    samples = pd.read_csv(options.samples, sep="\t", header=0, index_col=0)

    # Loading the barcodes
    #bcA = load_FASTA_barcode(options.bcA)
    bcB = load_FASTA_barcode(options.bcB)
    bcC = load_FASTA_barcode(options.bcC)
    bcD = load_FASTA_barcode(options.bcD)
    bcE = load_FASTA_barcode(options.bcE)

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

    # Producing the BC strains concordance scatter plots
    scatter_plot = pd.DataFrame()
    for sample_label in samples.index:
        # Counting strains per BC
        df = pd.read_csv(sample_label + "/" + sample_label + ".all_filter.tsv", sep="\t", index_col=0)
        df = df.groupby(["Barcode", "Strain"], as_index=False).UMI.count()
        df = df.pivot(index="Barcode", columns="Strain", values="UMI")
        df.loc[:, "Sample"] = sample_label
        scatter_plot = pd.concat([scatter_plot, df])
    scatter_plot.fillna(0, inplace=True)
    # Creating the FacetGrid
    g = sns.FacetGrid(scatter_plot, col="Sample", col_wrap=2, height=5)
    g = (g.map_dataframe(sns.scatterplot, "H1N1", "H3N2").set_titles("{col_name}"))
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle("Barcode strains concordance")
    plt.savefig("Barcode_strains_concordance.pdf", bbox_inches="tight")
    plt.clf()

    # Producing the reads/BC barplot
    barplot = pd.DataFrame()
    for sample_label in samples.index:
      df = pd.read_csv(sample_label + "/" + sample_label + ".droplet_total_reads.tsv", sep="\t", index_col=0)
      df = df.reset_index()
      df = df.groupby("Total_unique_reads_raw", as_index=False).count()
      df["Weights"] = df["Barcode"] * df["Total_unique_reads_raw"]
      df.loc[:, "Sample"] = sample_label
      barplot = pd.concat([barplot, df])
    # Creating the FacetGrid
    g = sns.FacetGrid(barplot, col="Sample", col_wrap=2, height=5)
    g = (g.map(weighted_hist, "Total_unique_reads_raw", "Weights", bins=np.logspace(0, 4, 100), kde=False).set_xlabels("#Reads/BC").set_ylabels("Count x #reads/BC").set(xscale="log").set_titles("{col_name}"))
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle("Read count per BC barplot")
    plt.savefig("Reads_per_BC_barplot.pdf", bbox_inches="tight")
    plt.clf()

    # Producing the strains distribution per RNA fragment stacked barplot
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    i = 0
    # Processing each sample data
    for sample_data in results:
        # Re-setting the index to be able to do a groupby
        df = sample_data[1].reset_index()
        df.set_index(["Strain", "RNA_fragment"], inplace=True)
        # Computing the percentages
        df = df.groupby("RNA_fragment").apply(lambda x: 100 * x / float(x.sum()))
        df.columns = ["Percentages"]
        # Resetting the index to be able to pivot
        df.reset_index(inplace=True)
        df = df.pivot(index="RNA_fragment", columns="Strain", values="Percentages")
        # Plotting
        df.plot.bar(ax=axes[sample_loc[sample_data[0]][0], sample_loc[sample_data[0]][1]], stacked=True, legend=False, title=sample_data[0], color=colors.values)
        i += 1
    # General figure management
    #fig.legend(handles=[mpatches.Patch(color=seg_colors["H1N1"], label="H1N1"), mpatches.Patch(color=seg_colors["H3N2"], label="H3N2"), mpatches.Patch(color=seg_colors["Undefined"], label="Undefined")])
    fig.legend(handles=leg)
    plt.tight_layout()
    plt.savefig("Strains_distribution_per_RNA_fragment.pdf", format="pdf")

    """
    # Clustermaps
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    i = 0
    # Processing each sample data
    for sample_label in samples.index:
        # Re-setting the index to be able to do a groupby
        df = pd.read_csv(sample_label + "/" + sample_label + ".clustermap.tsv", sep="\t")
        df.set_index(["Strain", "RNA_fragment"], inplace=True)
        # Computing the percentages
        df = df.groupby("RNA_fragment").apply(lambda x: 100 * x / float(x.sum()))
        df.columns = ["Percentages"]
        # Resetting the index to be able to pivot
        df.reset_index(inplace=True)
        df = df.pivot(index="RNA_fragment", columns="Strain", values="Percentages")
        # Plotting
        df.plot.bar(ax=axes[sample_loc[sample_data[0]][0], sample_loc[sample_data[0]][1]], stacked=True, legend=False, title=sample_data[0], color=colors.values)
        i += 1
    # General figure management
    # fig.legend(handles=[mpatches.Patch(color=seg_colors["H1N1"], label="H1N1"), mpatches.Patch(color=seg_colors["H3N2"], label="H3N2"), mpatches.Patch(color=seg_colors["Undefined"], label="Undefined")])
    fig.legend(handles=leg)
    plt.tight_layout()
    plt.savefig("Strains_distribution_per_RNA_fragment.pdf", format="pdf")


    clustermap = clustermap.loc[:, clustermap.sum(axis=0) != 0]  # Removing the columns that sum to 0 (FloatingPointError otherwise)
    # Normalizing the data by row and by col
    clustermap = clustermap.apply(lambda x: x/x.max(), axis=0)
    clustermap = clustermap.apply(lambda x: x/x.max(), axis=1)
    # Drawing and saving the clustermap
    sys.setrecursionlimit(10000)
    sns.set(font_scale=0.7)
    #ax = sns.clustermap(clustermap, col_cluster=False, vmax=25, cmap="magma_r", z_score=1)  # Drawing the clustermap
    ax = sns.clustermap(clustermap, col_cluster=False, cmap="magma_r")  # Drawing the clustermap
    plt.setp(ax.ax_heatmap.get_yticklabels(), rotation=0)  # Rotating the y labels
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=30)  # Rotating the x labels
    ax.ax_row_dendrogram.set_visible(False)  # Hiding the dendrogram
    plt.savefig(output_dir + sample_label + ".clustermap.pdf", format="pdf")
    plt.clf()
    sys.setrecursionlimit(1000)
    """

    # Producing the summary strains barplot
    #all_samples_mapping.groupby(["Strain"]).sum().T.plot(kind="bar", stacked=True, colormap="Set2", legend="reverse").legend(bbox_to_anchor=(1.5, 1))  # Not normalized
    all_samples_mapping.groupby(["Strain"]).sum().apply(lambda x: 100 * x / x.sum()).T.plot(kind="bar", stacked=True, colormap="Set2", legend="reverse").legend(bbox_to_anchor=(1.5, 1))
    plt.title("Strains distribution within samples")
    plt.savefig("Strains_barplot.pdf", bbox_inches="tight")
    plt.clf()
    # Producing the amplicons barplot
    all_samples_mapping.reset_index(inplace=True)
    amplicons_barplot = all_samples_mapping.loc[(all_samples_mapping["Strain"] != "Undefined") & (all_samples_mapping["Strain"] != "MappingAmbiguous") & (all_samples_mapping["RNA_fragment"] != "Undefined") & (all_samples_mapping["RNA_fragment"] != "MappingAmbiguous")].copy()
    amplicons_barplot = pd.melt(amplicons_barplot, id_vars=["Strain", "RNA_fragment"], value_vars=["i04", "i05", "i06", "i08"], var_name="Sample")
    amplicons_barplot["Amplicon"] = amplicons_barplot["Strain"] + "_" + amplicons_barplot["RNA_fragment"]
    ax = sns.barplot(data=amplicons_barplot, x="Amplicon", y="value", hue="Sample", ci=None)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title("Amplicons distribution within samples")
    plt.tight_layout()
    plt.savefig("Amplicons_barplot.pdf")
    plt.clf()
    # Producing the undefined amplicons barplot
    amplicons_barplot2 = all_samples_mapping.loc[(all_samples_mapping["Strain"] == "Undefined") | (all_samples_mapping["Strain"] == "MappingAmbiguous") | (all_samples_mapping["RNA_fragment"] == "Undefined") | (all_samples_mapping["RNA_fragment"] == "MappingAmbiguous")].copy()
    amplicons_barplot2 = pd.melt(amplicons_barplot2, id_vars=["Strain", "RNA_fragment"], value_vars=["i04", "i05", "i06", "i08"], var_name="Sample")
    amplicons_barplot2["Amplicon"] = amplicons_barplot2["Strain"] + "_" + amplicons_barplot2["RNA_fragment"]
    ax = sns.barplot(data=amplicons_barplot2, x="Amplicon", y="value", hue="Sample", ci=None)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title("Undefined amplicons distribution within samples")
    plt.tight_layout()
    plt.savefig("Amplicons_barplot.undefined.pdf")
    plt.clf()