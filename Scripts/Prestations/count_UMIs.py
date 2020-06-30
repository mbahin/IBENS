#!/usr/bin/python

"""
Copyright Mathieu Bahin, Auguste Genovesio and IBENS bioinformatics platform 2017

mathieu.bahin@biologie.ens.fr
auguste.genovesio@biologie.ens.fr

This software is a computer program whose purpose is to demultiplex
complex sequencing data.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""

__author__ = 'bahin'
''' Script to decipher how many different UMIs there were in each well of a sequencing plate. Written for DemultiplexHepCap_AB2016. '''

import argparse
import sys
import re
import editdistance
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import subprocess
import collections

# TODO
#   - UMI flex more generic?? re.match??
#   - Unbalanced output for plate barcodes (should be treated as col, row and UMI)
#   - Integrate the heatmaps??
#   - col/row_pattern useful?? (in config.txt)

### Functions

def compare_barcode(seq, barcode_start, barcode_stop, barcodes_list, allowed_dist=1, shift=0):
    """ Compare a part of a sequence (the barcode) to a list of expected barcodes. """
    # Defining the regular barcode part in the sequence
    barcode_seq = seq[(barcode_start + shift):(barcode_stop + shift)]
    if shift != 0:
        # If the (right) shift sets the barcode_stop further than the sequence end
        if (barcode_stop + shift) > len(seq):
            barcode_seq = seq[(barcode_start + shift):] # The start is shifted but the end is the end of the sequence
        # If the (left) shift sets the barcode_start before 0
        if (barcode_start + shift) < 0:
            barcode_seq = seq[0:(barcode_stop + shift)]
    # Comparing this barcode to the list of expected ones
    for barcode in barcodes_list:
        #barcode_exp = barcode[barcode_start:barcode_stop]
        barcode_exp = barcodes_list[barcode][barcode_start:barcode_stop]
        # Computing the distance between the barcode of the sequence and one of the list
        dist = editdistance.eval(barcode_exp, barcode_seq)
        if dist <= allowed_dist:
            # If the distance is lesser than the allowed one, then the barcode is assigned to the sequence
            return barcode
    # Otherwise the function returns "Unknown"
    return 'Unknown'


def multi_search_barcode(reads_file, barcode_start, barcode_stop, barcode_type, barcodes_list, allowed_dist, shift=0):
    """ Search for barcodes in a sequences in a specific position with more or less errors and shifts. """
    for read_ID, seq, qual in FastqGeneralIterator(open(reads_file)):
        # Getting the read ID
        read_ID = read_ID.split(' ')[0]
        if read_ID not in seqs:
            seqs[read_ID] = {'gene': 'noGene', 'col': 'Unknown', 'row': 'Unknown', 'UMI': 'Unknown', 'plate': 'Unknown', 'matching_plate_indexes': False}
        # The barcode is first checked at its regular location (from row_barcode_start to row_barcode_end)
        seqs[read_ID][barcode_type] = compare_barcode(seq, barcode_start, barcode_stop, barcodes_list)
        # And then, if a shift is allowed, shifted to the right and to the left
        if shift > 0:
            # If no match, try with a right shift
            if seqs[read_ID][barcode_type] == 'Unknown':
                seqs[read_ID][barcode_type] = compare_barcode(seq, barcode_start, barcode_stop, barcodes_list, allowed_dist, shift)
                # If no match, try with a left shift
                if seqs[read_ID][barcode_type] == 'Unknown':
                    seqs[read_ID][barcode_type] = compare_barcode(seq, barcode_start, barcode_stop, barcodes_list, allowed_dist, shift * -1)

### end Functions

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-index1', dest='index1', required=True, help="'Index1' file")
parser.add_argument('-index2', dest='index2', required=True, help="'Index2' file")
parser.add_argument('-read1', dest='read1', required=True, help="Trimmed 'Read1' file")
parser.add_argument('-c', dest='config', required=True, help='Plate column and row indexes')
parser.add_argument('-o', dest='output', required=True, help='Output files basename')
parser.add_argument('-f', dest='failed_basename', help='Failed files basename')
options = parser.parse_args()

# Indexing the plate information
print "##### Running the script #####"
print "Indexing the plate column and row indexes..."
genes = collections.OrderedDict()
wells = collections.OrderedDict()  # {<plate_ID>: {<col_ID>: {<row_ID>: {<UMI_seq>: {<gene>: count}}}}
cols_rows_sets = collections.OrderedDict()  # {<set_ID>: ["Rxx/Cxx", ..., "Rxx/Cxx"]}
plates = collections.OrderedDict()  # {<plate_ID>: {"seq": XXX, "matching": ["Cxx", ..., "Cxx", "Rxx", ..., "Rxx"]}}
cols = collections.OrderedDict()  # {<col_ID>: <col_seq>}
rows = collections.OrderedDict()  # {<row_ID>: <row_seq>}
config = {}
tag = "general"  # Marker to follow the progress in the config file (general > gene > col > row)
with open(options.config, "r") as config_file:
    for line in config_file:
        if (tag == "general") and (line.rstrip() != "##Genes"):  # Getting the general config variables
            config[line.split("\t")[0]] = line.rstrip().split("\t")[1]
        elif (tag == "general") and (line.rstrip() == "##Genes"):  # First line, genes part
            tag = "gene"
        elif (tag == "gene") and (line.rstrip() != "##Cols/Rows sets"):  # Indexing the genes
            genes[line.split("\t")[0]] = line.rstrip().split("\t")[1]
        elif (tag == "gene") and (line.rstrip() == "##Cols/Rows sets"):  # Getting to the cols/rows sets part
            tag = "sets"
        elif (tag == "sets") and (line.rstrip() != "##Plates"):  # Indexing the expected cols/rows sets
            cols_rows_sets[line.split("\t")[0]] = line.strip().split("\t")[1].split(",")
        elif (tag == "sets") and (line.rstrip() == "##Plates"):  # Getting to the plates part
            tag = "plate"
        elif (tag == "plate") and (line.rstrip() != "##Columns"):  # Indexing the plates
            cols_rows_sets_IDs = []
            sets_list = [ids.split(",") for ids in line.rstrip().split("\t")[2:]]  # List of the matching col and row barcodes per plate barcode
            for set_list in sets_list:
                new_set = []
                for col_row_id in set_list:  # Transforming col/row set ID into col and row IDs list
                    new_set.append(cols_rows_sets[col_row_id])
                cols_rows_sets_IDs.append(new_set)
            plates[line.split("\t")[0]] = {"seq": line.split("\t")[1], "matching": cols_rows_sets_IDs}  # Registering the expected row and col IDs for the plate
        elif (tag == "plate") and (line.rstrip() == "##Columns"):  # Getting to the columns part
            tag = "col"
        elif (tag == "col") and (line.rstrip() != "##Rows"):  # Indexing the columns
            cols[line.split("\t")[0]] = line.rstrip().split("\t")[1]
        elif (tag == "col") and (line.rstrip() == "##Rows"):  # Getting to the rows part
            tag = "row"
        elif tag == "row":  # Indexing the rows
            rows[line.split("\t")[0]] = line.rstrip().split("\t")[1]

# Checking the distances between col barcodes and between row barcodes
for feat in ((cols, 'column'), (rows, 'row')):
    min_dist = -1
    for ID in feat[0]:  # First instance of the barcodes
        for ID2 in feat[0]:  # Second instance of the barcodes
            if ID != ID2:  # We don't want to compare the same barcode to himself
                if feat[0][ID] == feat[0][ID2]:
                    sys.exit('Error: there are 2 similar barcodes (seq: ' + feat[0][ID] + ') in the ' + feat[1] + ' barcodes list. Aborting.')
                else:
                    dist = int(editdistance.eval(feat[0][ID], feat[0][ID2]))
                # Updating the minimum distance
                if (min_dist == -1) or (dist < min_dist):
                    min_dist = dist
                    if min_dist <= (int(config['indexes_barcode_error']) * 2):
                        sys.exit('Error: there is at least one pair of ' + feat[1] + ' barcodes that have a distance of ' + str(min_dist) + ' although ' + config['indexes_barcode_error'] + ' error(s) is/are allowed on the barcode. This is inconsistent. Aborting.')

# Initializing the plate with only the expected col and row IDs for each plate ID
for plate_ID in plates:
    wells[plate_ID] = collections.OrderedDict()
    for col_ID in list(set([c for double_sublist in plates[plate_ID]["matching"] for sublist in double_sublist for c in sublist]).intersection(cols)):
        wells[plate_ID][col_ID] = collections.OrderedDict()
        for row_ID in list(set([r for double_sublist in plates[plate_ID]["matching"] for sublist in double_sublist for r in sublist]).intersection(rows)):
            wells[plate_ID][col_ID][row_ID] = collections.OrderedDict()

# Managing the pattern boundaries
for pat in re.finditer('XXXXXX', config['col_pattern']):
    col_pat_beg = int(pat.start())
    col_pat_end = int(pat.end())
for pat in re.finditer('XXXXXX', config['row_pattern']):
    row_pat_beg = int(pat.start())
    row_pat_end = int(pat.end())

if ((col_pat_beg - int(config['indexes_shift'])) + int(config['indexes_barcode_error'])) < 0:
    sys.exit('Error: The column barcode is too close to the start of the sequence to apply the indexes shift. Aborting.')
elif (col_pat_beg - int(config['indexes_shift'])) < 0:
    sys.stderr.write('Warning: the column barcode is too close to the start of the sequence to properly apply the indexes shift, systematic error(s) will be taken into account instead of the overhanging bp.\n')
if ((row_pat_beg - int(config['indexes_shift'])) + int(config['indexes_barcode_error'])) < 0:
    sys.exit('Error: The row barcode is too close to the start of the sequence to apply the indexes shift. Aborting.')
elif (row_pat_beg - int(config['indexes_shift'])) < 0:
    sys.stderr.write('Warning: the row barcode is too close to the start of the sequence to properly apply the shift, systematic error(s) will be taken into account instead of the overhanging bp.\n')
if ((col_pat_end + int(config['indexes_shift'])) - int(config['indexes_barcode_error'])) > len(config['col_pattern']):
    sys.exit('Error: The column barcode is too close to the end of the sequence to properly apply the indexes shift. Aborting.')
elif (col_pat_end + int(config['indexes_shift'])) > len(config['col_pattern']):
    sys.stderr.write('Warning: the column barcode is too close to the end of the sequence to properly apply the shift, systematic error(s) will be taken into account instead of the overhanging bp.\n')
if ((row_pat_end + int(config['indexes_shift'])) - int(config['indexes_barcode_error'])) > len(config['row_pattern']):
    sys.exit('Error: The row barcode is too close to the end of the sequence to apply the indexes shift. Aborting.')
elif (row_pat_end + int(config['indexes_shift'])) > len(config['row_pattern']):
    sys.stderr.write('Warning: the row barcode is too close to the end of the sequence to properly apply the indexes shift, systematic error(s) will be taken into account instead of the overhanging bp.\n')

# Processing 'Index1' file (columns)
seqs = {}  # {<read_ID>: {'gene': tooShortAfter3TrimUv/noGene/<gene_name>, 'col': <col_ID>/Unknown, 'row': <row_ID>/Unknown, UMI: <UMI_seq>/Unknown, 'plate': <plate_ID>/Unknown, matching_plate_indexes: True/False}}
print 'Processing the column indexes file...'
multi_search_barcode(options.index1, col_pat_beg, col_pat_end, 'col', cols, int(config['indexes_barcode_error']), int(config['indexes_shift']))

# Processing 'Index2' file (rows)
print 'Processing the row indexes file...'
multi_search_barcode(options.index2, row_pat_beg, row_pat_end, 'row', rows, int(config['indexes_barcode_error']), int(config['indexes_shift']))

# Processing 'Read1' file (UMI + plate)
print 'Processing the "Read1" file...'

# Trimming the Uv sequence on the 3' end of the sequences (keep only sequences at least 185nt long)
with open(os.devnull, 'w') as devnull, open('read1.3TrimUv.fastq', 'w') as output_file:
    p = subprocess.Popen(['cutadapt', '-a', config['Uv'], '-m', config['seq_min_length'], '--untrimmed-output', options.failed_basename + 'untrimmedUv.fasta', '--too-short-output', options.failed_basename + 'tooShortAfter3TrimUv.fasta', '-e', config['Uv_error'], options.read1], stdout=output_file, stderr=devnull, close_fds=True)
    p.wait()  # Wait for the cutadapt command to be finished
# Getting the minimum distances between genes within the list
min_gene_dist = len(genes[min(genes, key=genes.get)])  # Initialization of the maximum error allowed on the gene with the minimum gene length
for g1 in genes:
    for g2 in genes:
        if g1 != g2:
            min_gene_dist = min(min_gene_dist, editdistance.eval(genes[g1], genes[g2]))
min_genes_dist = int(min_gene_dist / 2) - (min_gene_dist % 2 == 0)  # Getting the maximum distance allowed: one less than half (if even number, floor rounded half otherwise)
# Catching the UMI (with one mistake allowed)
for read_ID, seq, qual in FastqGeneralIterator(open('read1.3TrimUv.fastq')):
    # Getting the read ID
    read_ID = read_ID.split(' ')[0]
    # Matching the UMI pattern (exactly and then with one error allowed)
    for p in (r'...AG..AG...AG..', r'....G..AG...AG..', r'...A...AG...AG..', r'...AG...G...AG..', r'...AG..A....AG..', r'...AG..AG....G..', r'...AG..AG...A...'):
        match = re.match(p, seq[-16:])
        if match:
            seqs[read_ID]['UMI'] = '-'.join((match.group(0)[0:3], match.group(0)[5:7], match.group(0)[9:12], match.group(0)[14:16]))
            break
    #seqs[read_ID]['UMI'] = "NoUMI"
    # Identifying the gene
    for g in genes:
        g_seq = seq[4:len(genes[g])]
        if editdistance.eval(g_seq, genes[g][4:len(genes[g])]) <= min_genes_dist:
            seqs[read_ID]['gene'] = g
            break
# Remove the temporary trimmed file
os.remove('read1.3TrimUv.fastq')

# Indexing the untrimmed reads
report = {'untrimmedUv': []}  # Container for all the failed reads
for seq_record in SeqIO.parse(options.failed_basename + 'untrimmedUv.fasta', "fasta"):
    report['untrimmedUv'].append(seq_record.id)
    seqs[seq_record.id]['gene'] = 'untrimmedUv'

# Indexing the too short sequences after the Uv sequence trimming
report['tooShortAfter3TrimUv'] = []
for seq_record in SeqIO.parse(options.failed_basename + 'tooShortAfter3TrimUv.fasta', "fasta"):
    report['tooShortAfter3TrimUv'].append(seq_record.id)
    seqs[seq_record.id]['gene'] = 'tooShortAfter3TrimUv'

# 5' trimming the Uv sequence on read1 to get the plate barcode
with open(os.devnull, 'w') as devnull, open('read1.5TrimUv.fastq', 'w') as output_file:
    p = subprocess.Popen(['cutadapt', '-g', config['Uv'], '--discard-untrimmed', '-e', config['Uv_error'], options.read1], stdout=output_file, stderr=devnull, close_fds=True)
    p.wait()  # Wait for the cutadapt command to be finished
# Catching the plate barcode (with one mistake allowed)
plates_ID_seq = {}
for plate_ID in plates:
    plates_ID_seq[plate_ID] = plates[plate_ID]['seq']
multi_search_barcode('read1.5TrimUv.fastq', 3, 9, 'plate', plates_ID_seq, int(config['plate_barcode_error']), int(config['plate_shift']))
# Remove the temporary trimmed file
os.remove('read1.5TrimUv.fastq')

# Reordering the data
for cat in ('exploited', 'trashed', 't_plate', 't_conflict', 't_UMI', 't_col', 't_row', 't_UMI_col', 't_UMI_row', 't_col_row', 't_UMI_col_row'):
    report[cat] = {}
    for gene in genes:
        report[cat][gene] = []
total = 0
noGene_count = []
gene_count = []

for read_ID in seqs:
    total += 1
    if seqs[read_ID]['gene'] in ['tooShortAfter3TrimUv', 'untrimmedUv']:
        continue
    if seqs[read_ID]['gene'] != 'noGene':
        gene_count.append(read_ID)
        gene_name = seqs[read_ID]['gene']
        if seqs[read_ID]['plate'] == 'Unknown':
            report['t_plate'][gene_name].append(read_ID)  # No plate
            report['trashed'][gene_name].append(read_ID)
        else:
            # Checking if the plate barcode is matching the indexes it is supposed to (by default it's False, we check if it's True)
            if (seqs[read_ID]["col"] != "Unknown") and (seqs[read_ID]["row"] != "Unknown"):
                for cols_rows_set in plates[seqs[read_ID]["plate"]]["matching"]:
                    if (seqs[read_ID]["col"] in [x for sublist in cols_rows_set for x in sublist]) and (seqs[read_ID]["row"] in [x for sublist in cols_rows_set for x in sublist]):
                        seqs[read_ID]["matching_plate_indexes"] = True
                        break
                else:
                    seqs[read_ID]["matching_plate_indexes"] = False
            elif seqs[read_ID]["col"] != "Unknown":
                for cols_rows_set in plates[seqs[read_ID]["plate"]]["matching"]:
                    if seqs[read_ID]["col"] in [x for sublist in cols_rows_set for x in sublist]:
                        seqs[read_ID]["matching_plate_indexes"] = True
                        break
                else:
                    seqs[read_ID]['matching_plate_indexes'] = False
            elif seqs[read_ID]["row"] != "Unknown":
                for cols_rows_set in plates[seqs[read_ID]["plate"]]["matching"]:
                    if seqs[read_ID]["row"] in [x for sublist in cols_rows_set for x in sublist]:
                        seqs[read_ID]["matching_plate_indexes"] = True
                        break
                else:
                    seqs[read_ID]["matching_plate_indexes"] = False
            else:
                seqs[read_ID]["matching_plate_indexes"] = True
            if not seqs[read_ID]['matching_plate_indexes']:
                report['t_conflict'][gene_name].append(read_ID)  # Plate barcode not matching indexes barcodes
                report['trashed'][gene_name].append(read_ID)
            else:
                if (seqs[read_ID]['col'] != 'Unknown') and (seqs[read_ID]['row'] != 'Unknown') and (seqs[read_ID]['UMI'] != 'Unknown'):
                    report['exploited'][gene_name].append(read_ID)
                    if seqs[read_ID]['UMI'] not in wells[seqs[read_ID]['plate']][seqs[read_ID]['col']][seqs[read_ID]['row']]:
                        # Initialization of the genes dicts for this UMI
                        wells[seqs[read_ID]['plate']][seqs[read_ID]['col']][seqs[read_ID]['row']][seqs[read_ID]['UMI']] = {}
                        for gene in genes:
                            wells[seqs[read_ID]['plate']][seqs[read_ID]['col']][seqs[read_ID]['row']][seqs[read_ID]['UMI']][gene] = []
                    # Registering the read
                    wells[seqs[read_ID]['plate']][seqs[read_ID]['col']][seqs[read_ID]['row']][seqs[read_ID]['UMI']][gene_name].append(read_ID)
                else:
                    if (seqs[read_ID]['UMI'] == 'Unknown') and (seqs[read_ID]['col'] != 'Unknown') and (seqs[read_ID]['row'] != 'Unknown'):  # No UMI
                        report['t_UMI'][gene_name].append(read_ID)
                    if (seqs[read_ID]['col'] == 'Unknown') and (seqs[read_ID]['row'] != 'Unknown') and (seqs[read_ID]['UMI'] != 'Unknown'):  # No column
                        report['t_col'][gene_name].append(read_ID)
                    if (seqs[read_ID]['row'] == 'Unknown') and (seqs[read_ID]['col'] != 'Unknown') and (seqs[read_ID]['UMI'] != 'Unknown'):  # No row
                        report['t_row'][gene_name].append(read_ID)
                    if (seqs[read_ID]['UMI'] == 'Unknown') and (seqs[read_ID]['col'] == 'Unknown') and (seqs[read_ID]['row'] != 'Unknown'):  # No UMI and column
                        report['t_UMI_col'][gene_name].append(read_ID)
                    if (seqs[read_ID]['UMI'] == 'Unknown') and (seqs[read_ID]['row'] == 'Unknown') and (seqs[read_ID]['col'] != 'Unknown'):  # No UMI and row
                        report['t_UMI_row'][gene_name].append(read_ID)
                    if (seqs[read_ID]['col'] == 'Unknown') and (seqs[read_ID]['row'] == 'Unknown') and (seqs[read_ID]['UMI'] != 'Unknown'):  # No column and row
                        report['t_col_row'][gene_name].append(read_ID)
                    if (seqs[read_ID]['UMI'] == 'Unknown') and (seqs[read_ID]['col'] == 'Unknown') and (seqs[read_ID]['row'] == 'Unknown'):  # No UMI, column and row
                        report['t_UMI_col_row'][gene_name].append(read_ID)

                    report['trashed'][gene_name].append(read_ID)
    else:
        noGene_count.append(read_ID)

# Printing the output files
print 'Printing the output files...'
# Writing the report file
with open(options.output + '.report.tsv', 'w') as report_file:
    report_file.write('##### Report #####\n')
    report_file.write('Total sequences\t' + str(total) + '\n')
    report_file.write('Untrimmed reads (Uv not found)\t' + str(len(report['untrimmedUv'])) + '\n')
    report_file.write('Rotten reads (too short after trimming the Uv)\t' + str(len(report['tooShortAfter3TrimUv'])) + '\n')
    report_file.write('Gene not matched sequences\t' + str(len(noGene_count)) + '\n')
    report_file.write('Gene matched sequences\t' + str(len(gene_count)) + '\n\n')
    for gene in genes:
        report_file.write('\t' + gene)
    report_file.write('\tTotal\n')
    # Exploited/unexploited sequences
    for cat in (('Exploited sequences', report['exploited']), ('Unexploited sequences', report['trashed'])):
        report_file.write(cat[0])
        for gene in genes:
            report_file.write('\t' + str(len(cat[1][gene])))
        report_file.write('\t' + str(sum(len(lst) for lst in cat[1].values())) + '\n')
    # Unexploited sequences details
    for cat in (('plate barcode', report['t_plate']), ('plate/indexes barcodes conflict', report['t_conflict']), ('UMI sequence', report['t_UMI']), ('Index1 (col) barcode', report['t_col']), ('Index2 (row) barcode', report['t_row']), ('UMI sequence and Index1 barcode', report['t_UMI_col']), ('UMI and Index2 barcode', report['t_UMI_row']), ('Index1 and Index2 barcodes', report['t_col_row']), ('UMI sequence, Index1 and Index2 barcodes', report['t_UMI_col_row'])):
        report_file.write('  Discarded because of the ' + cat[0])
        for gene in genes:
            report_file.write('\t' + str(len(cat[1][gene])))
        report_file.write('\t' + str(sum(len(lst) for lst in cat[1].values())) + '\n')
# Writing the summarized output file
with open(options.output + '.tsv', 'w') as output_file:
    output_file.write('Plate\tPlate barcode\tColumn\tColumn barcode\tRow\tRow barcode')
    for gene in genes:
        output_file.write('\tUMI count in ' + gene + '\tAverage count per UMI in ' + gene)
    output_file.write('\n')
    for plate_ID in plates:
        for cols_rows_set in plates[plate_ID]["matching"]:
            for col_ID in cols_rows_set[0]:
                for row_ID in cols_rows_set[1]:
                    output_file.write('\t'.join((plate_ID, str(plates[plate_ID]['seq']), col_ID, str(cols[col_ID]), row_ID, str(rows[row_ID]))))
                    for gene in genes:
                        nb = 0
                        counts = 0
                        for UMI in wells[plate_ID][col_ID][row_ID]:
                            if len(wells[plate_ID][col_ID][row_ID][UMI][gene]) > 0:
                                nb += 1
                                counts += len(wells[plate_ID][col_ID][row_ID][UMI][gene])
                        if nb > 0:
                            output_file.write('\t' + str(nb) + '\t' + str(round(counts / nb, 2)))
                        else:
                            output_file.write('\t0\tNA')
                    output_file.write('\n')
# Writing the detailed output file
with open(options.output + '.detailed.tsv', 'w') as output_detailed_file:
    output_detailed_file.write('Plate\tPlate barcode\tColumn\tColumn barcode\tRow\tRow barcode\tUMI')
    for gene in genes:
        output_detailed_file.write('\tCount in ' + gene)
    output_detailed_file.write('\n')
    for plate_ID in plates:
        for cols_rows_set in plates[plate_ID]["matching"]:
            for col_ID in cols_rows_set[0]:
                for row_ID in cols_rows_set[1]:
                    for UMI in wells[plate_ID][col_ID][row_ID]:
                        output_detailed_file.write('\t'.join((plate_ID, str(plates[plate_ID]['seq']), col_ID, str(cols[col_ID]), row_ID, str(rows[row_ID]), UMI)))
                        for gene in genes:
                            output_detailed_file.write('\t' + str(len(wells[plate_ID][col_ID][row_ID][UMI][gene])))
                        output_detailed_file.write('\n')

# Reporting the reads that failed if asked
if options.failed_basename:
    print "Printing the failed report files..."
    # Reads that failed because no gene was found
    with open(options.failed_basename + "Gene.txt", "w") as failed_file:
        for read_ID in noGene_count:
            failed_file.write(read_ID + "\n")
    # A gene was found but failed for another reason (UMI, col or row combination)
    for cat in (("Plate", report["t_plate"]), ("Barcodes_conflict", report["t_conflict"]), ("UMI", report["t_UMI"]), ("Index1", report["t_col"]), ("Index2", report["t_row"]), ("UMI_Index1", report["t_UMI_col"]), ("UMI_Index2", report["t_UMI_row"]), ("Index1_Index2", report["t_col_row"]), ("UMI_Index1_Index2", report["t_UMI_col_row"])):
        with open(options.failed_basename + cat[0] + ".txt", "w") as failed_file:
            for gene in cat[1]:
                for read_ID in cat[1][gene]:
                    failed_file.write(read_ID + "\n")
print "##### End of script #####"