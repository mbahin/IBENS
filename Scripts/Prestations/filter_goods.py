#!/usr/bin/python

__author__ = 'bahin'
''' Script to keep only the reads that have a proper column/row index pair matching with the plate barcodes. '''

import argparse
import sys
from Bio import SeqIO
import editdistance

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', required=True, help='Fastq input file.')
parser.add_argument('-o', dest='output', required=True, help='Filtered fastq output file.')
options = parser.parse_args()

# Indexing the plate barcodes
plate_barcodes = {}
plate_barcodes['P41'] = ['R01', 'R02', 'R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26']
plate_barcodes['P44'] = ['R01', 'R02', 'R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38']
plate_barcodes['P47'] = ['R09', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26']
plate_barcodes['P43'] = ['R09', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38']

# Parsing the input file and writing the output one
with open(options.input, 'r') as input_file, open(options.output, 'w') as output_file:
    for line in input_file:
        pbc = line.split(' ')[2]
        if (pbc != 'Unknown') and (line.split(' ')[3] in plate_barcodes[pbc]) and (line.split(' ')[4] in plate_barcodes[pbc]):
            output_file.write(line.rstrip() + '\n')