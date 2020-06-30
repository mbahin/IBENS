#!/usr/bin/python

__author__ = "bahin"
""" Script to check whether there really are way more substitions in the barcodes than indels in DemultiplexInfluenza data. """

import argparse
from Bio import SeqIO
import sys

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument("--sample", dest="sample", required=True, help="Reads file")
parser.add_argument("-a", dest="bcA", default="/projects/biocompan/Bioinformatics_services/201810_Pasteur_Isel-Griffith_DemultiplexInfluenza_A2018/Data/A_forward.fasta", help="BarcodeA FASTA file")
options = parser.parse_args()

# Initializations
alphabet = ["A", "T", "G", "C"]
exact = {}
hamming = {}
levenshtein = {}
other = {}

# Browsing barcodes file and creating the Hamming and Levenshtein dicts
for record in SeqIO.parse(options.bcA, "fasta"):
    bc = record.seq
    exact[str(bc)] = 0
    # Hamming distance
    for pos in range(len(bc)):
        for nt in alphabet:
            if bc[pos] != nt:
                #print "Sub", bc[0:pos] + nt + bc[pos + 1:]
                hamming[bc[0:pos] + nt + bc[pos + 1:]] = 0

    # Levenshtein distance
    for pos in range(len(bc)):
        for nt in alphabet:
            # Deletion
            levenshtein[bc[0:pos] + bc[pos + 1:] + nt] = 0
            #print "Del", bc[0:pos] + bc[pos + 1:] + nt
            # Insertion
            levenshtein[bc[0:pos] + nt + bc[pos:-1]] = 0
            #print "In", bc[0:pos] + nt + bc[pos:-1]
# Browsing the data file to identify the barcodes
exact_cpt = 0
hamming_cpt = 0
levenshtein_cpt = 0
other_cpt = 0
#i = 0
for record in SeqIO.parse(options.sample, "fastq"):
    #i += 1
    #if i % 10000 == 0:
        #print str(i) + " processed reads..."
    barcode = record.seq[0:16]
    if barcode in exact:
        exact_cpt += 1
        exact[barcode] += 1
    elif barcode in hamming:
        hamming_cpt += 1
        hamming[barcode] += 1
    elif barcode in levenshtein:
        levenshtein_cpt += 1
        levenshtein[barcode] += 1
    else:
        other_cpt += 1
        if barcode in other:
            other[barcode] += 1
        else:
            other[barcode] = 1

# Display results
print str(exact_cpt) + " exact matches."
print str(hamming_cpt) + " matches with 1 Hamming distance."
print str(levenshtein_cpt) + " matches with 1 Levenshtein distance."
print "Other matches:", other_cpt
print "Other matches:", other