#!/usr/bin/python

__author__ = 'bahin'
''' Script to filter a fastq file. The input file can be filtered either (exclusively) on an IDs list or on a sequence (with a distance). If the file is filtered on a sequence, the script can output the list of kept IDs. '''

import argparse
import sys
from Bio import SeqIO
import editdistance

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', required=True, help='Fastq input file.')
parser.add_argument('-r', dest='remove', action='store_true', help='Flag to remove the matching reads (with sequence or ID) instead of keeping them (default behaviour).')
parser.add_argument('-s', dest='seq', help='Sequence to look for in the input file.')
parser.add_argument('-d', dest='dist', type=int, help='Distance threshold between the sequence to look for and the reads for them to be kept.')
parser.add_argument('-l', dest='list', help='IDs list on which to filter the input file.')
parser.add_argument('-oI', dest='output_IDs', help='Output file with the kept IDs (if the input file was filtered on a sequence).')
parser.add_argument('-o', dest='output', required=True, help='Filtered fastq output file.')
options = parser.parse_args()

# Checking arguments
if not (bool(options.seq) ^ bool(options.list)):
    sys.exit('One and only one of the arguments "-s" and "-l" must be set. Aborting.')
if options.list and options.dist:
    sys.exit('Arguments "-l" and "-d" were set. Irrelevant. Aborting.')
if options.list and options.output_IDs:
    sys.exit('Arguments "-l" and "-oI" were set. Irrelevant. Aborting.')

# If the list mode is used, indexing the list
if options.list:
    input_IDs = []
    with open(options.list, 'r') as IDs_file:
        for line in IDs_file:
            input_IDs.append(line.rstrip())

# Parsing the input file
print 'Filtering the input file...'
i = 0
matching_reads = []
not_matching_reads = []
for seq_record in SeqIO.parse(options.input, "fastq"):
    i += 1
    if i % 1000000 == 0:
        sys.stdout.write(str(i) + ' sequences processed...\n')
        sys.stdout.flush()
    if (options.seq and (editdistance.eval(seq_record.seq, options.seq) <= options.dist)) or (options.list and (seq_record.id in input_IDs)):
        matching_reads.append(seq_record)
    else:
        not_matching_reads.append(seq_record)

# Writing the fastq output file
print 'Writing the output file(s)...'
if options.remove: # Remove matching reads
    to_write = not_matching_reads
else:
    to_write = matching_reads
SeqIO.write(to_write, options.output, "fastq")

# Writing the output IDs list (if needed)
if options.output_IDs:
    with open(options.output_IDs, 'w') as output_IDs_file:
        for seq_record in to_write:
            output_IDs_file.write(seq_record.id + '\n')