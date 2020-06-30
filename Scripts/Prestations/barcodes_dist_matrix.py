#!/usr/bin/python

__author__ = 'bahin'
''' Script to find the couples of sequences the most and least close (compared with Levenshtein distance) to each other within a set of sequences. '''

import argparse
import editdistance

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', required=True, help='<description>')
parser.add_argument('-o', dest='output', required=True, help='<description>')
options = parser.parse_args()

# Collecting the sequences from input file
seqs = []
with open(options.input, 'r') as input_file:
    for line in input_file:
        seqs.append(line.rstrip().upper())

# Computing and writing the distance matrix
dists = {}
with open(options.output, 'w') as output_file:
    output_file.write('\t' + '\t'.join(seqs) + '\n')
    for s1 in seqs:
        output_file.write(s1)
        for s2 in seqs:
            dist = int(editdistance.eval(s1, s2))
            if s1 != s2:
                if dist not in dists:
                    dists[dist] = set()
                dists[dist].add((min(s1, s2), max(s1, s2)))
            output_file.write('\t' + str(dist))
        output_file.write('\n')
print 'Min:', min(dists), dists[min(dists)]
print 'Max:', max(dists), dists[max(dists)]