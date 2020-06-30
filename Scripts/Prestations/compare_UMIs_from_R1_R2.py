#!/usr/bin/python

__author__ = 'bahin'
''' <description> '''

import argparse
import os
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-read1', dest='read1', required=True, help='<description>')
parser.add_argument('-read2', dest='read2', help='<description>')
parser.add_argument('-o', dest='output', required=True, help='<description>')
options = parser.parse_args()

# Preparing the results dict

# Getting UMI from R1
print 'Getting info from read1...'
res = {}
with open(os.devnull, 'w') as devnull, open('read1.3-trimmed.fastq', 'w') as output_file:
    p = subprocess.Popen(['cutadapt', '-a', 'GTGTAGGTCGTTCGCTCCAA', '-m', '185', '-e', '0.25', options.read1], stdout=output_file, stderr=devnull, close_fds=True)
    p.wait() # Wait for the cutadapt command to be finished
for read_id, seq, qual in FastqGeneralIterator(open('read1.3-trimmed.fastq')):
    # Getting the read IDs
    read_id = read_id.split(' ')[0]
    res[read_id] = {'R1': 'Unknown', 'R2': 'Unknown'}
    # Matching the UMI pattern (exactly and then with one error allowed)
    for p in (r'...AG..AG...AG..', r'....G..AG...AG..', r'...A...AG...AG..', r'...AG...G...AG..', r'...AG..A....AG..', r'...AG..AG....G..', r'...AG..AG...A...'):
        match = re.match(p, seq[-16:])
        if match:
            res[read_id]['R1'] = '-'.join((match.group(0)[0:3], match.group(0)[5:7], match.group(0)[9:12], match.group(0)[14:16]))
            break

# Getting UMI from R2
print 'Getting info from read2...'
for read_id, seq, qual in FastqGeneralIterator(open(options.read2)):
    # Getting the read IDs
    read_id = read_id.split(' ')[0]
    if read_id not in res:
        res[read_id] = {'R1': 'trimmed', 'R2': 'Unknown'}
    # Matching the UMI pattern (exactly and then with one error allowed)
    for p in (r'..CT...CT..CT...', r'...T...CT..CT...', r'..C....CT..CT...', r'..CT....T..CT...', r'..CT...C...CT...', r'..CT...CT...T...', r'..CT...CT..C....'):
        match = re.match(p, seq[0:16])
        if match:
            pat = '-'.join((match.group(0)[0:2], match.group(0)[4:7], match.group(0)[9:11], match.group(0)[13:16]))
            my_seq = Seq(pat, generic_dna)
            res[read_id]['R2'] = str(my_seq.reverse_complement())
            break
    ''' Version no mismatch on the triple 'CT'
    match = re.match(r'..CT...CT..CT...', seq[0:16])
    if match:
        pat = '-'.join((match.group(0)[0:2], match.group(0)[4:7], match.group(0)[9:11], match.group(0)[13:16]))
        my_seq = Seq(pat, generic_dna)
        res[read_id]['R2'] = str(my_seq.reverse_complement())
    '''

# Results report
with open(options.output, 'w') as output_file:
    output_file.write('Read ID\tFrom read1\tFrom read2\n')
    for read_id in res:
        output_file.write(read_id + '\t' + res[read_id]['R1'] + '\t' + res[read_id]['R2'] + '\n')