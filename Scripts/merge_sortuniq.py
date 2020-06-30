__author__ = 'bahin'
''' Script to create a table out of several files obtained with "sort uniq".'''

import argparse
import sys

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', nargs='+', required=True) # Input 'sort uniq' files
parser.add_argument('-o', dest='output', required=True)
options = parser.parse_args()

# Importing the 'common_func' module
sys.path.insert(1, '/data/biocomp/bahin/biocomp-mb-scripts/Scripts')
from common_func import init_dict
sys.path.remove('/data/biocomp/bahin/biocomp-mb-scripts/Scripts')

# Indexing 'sort uniq' files
stats = {}
max_occ = 0
for f in options.input:
    with open(f, 'r') as sortuniq_file:
        for line in sortuniq_file:
            init_dict(stats, line.rstrip().split('\t')[1], {})
            stats[line.rstrip().split('\t')[1]][options.input.index(f)] = line.split('\t')[0]
            max_occ = max(max_occ, int(line.rstrip().split('\t')[1]))

# Writing the output
with open(options.output, 'w') as output_file:
    output_file.write('\t' + '\t'.join((options.input)) + '\n')
    for i in range (1, max_occ + 1):
        output_file.write(str(i))
        for f in options.input:
            if (str(i) in stats) and (options.input.index(f) in stats[str(i)]):
                output_file.write('\t' + stats[str(i)][options.input.index(f)])
            else:
                output_file.write('\t0')
        output_file.write('\n')