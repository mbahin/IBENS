#!/usr/bin/python

__author__ = 'bahin'
''' Script to list all the DNA patterns to a distance of 0 to 1 to an input pattern. '''

import argparse

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='pattern', required=True, help='Input DNA pattern')
parser.add_argument('-o', dest='output', required=True, help='Output file with one sequence per line.')
options = parser.parse_args()

with open(options.output, 'w') as output_file:
    output_file.write(options.pattern + '\n') # First line is the pattern itself
    for pos in range(len(options.pattern)):
        beg = options.pattern[0:pos] # Part before the snp
        for n in ['A', 'T', 'G', 'C']: # Trying the 4 nt letters
            if n != options.pattern[pos]: # We skip the one already in the input pattern
                output_file.write(beg + n + options.pattern[pos + 1:] + '\n')
                seq = ''