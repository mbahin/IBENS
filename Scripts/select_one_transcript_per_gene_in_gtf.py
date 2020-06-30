__author__ = 'bastiane'

#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import re
import sys
import argparse
# =====================================================================================================================
# STRUCTURE D'UN FICHIER GTF :
# 0:chromosome    1:source    2:feature   3:start  4:end    5:score     6:strand    7:frame  8:attributes
# ---------------------------------------------------------------------------------------------------------------------
# Selects the longest transcript.

def best_tr(gtf_file):

    genes = {}
    cpt = 0

    gtf = open(gtf_file,'r')
    for line in gtf:
        split_line = re.split("\s+", line)
        if split_line[2] == "exon":
            last_elt = ("").join(split_line[8:])
            attributes = re.split('\s*;\s*',last_elt)

            gene = 0
            for attribute in attributes:
                gene = re.search('^gene_id"(\w+\.*\w*)"',attribute)
                if gene:
                    gene_id = gene.groups()[0]
                    if gene_id not in genes:
                        cpt += 1
                        if cpt%10000 == 0:
                            print(str(cpt)+" genes scanned")
                        genes[gene_id] = ("X", 0)

                    found_tr=""
                    found_nb=""
                    for attribute in attributes:
                        tr = re.search('^transcript_id"(\w+\.*\w*)"',attribute)
                        if tr:
                            found_tr=tr.groups()[0]
                        exon_nb = re.search('^exon_number"(\d+)"',attribute)
                        if exon_nb:
                            found_nb=int(exon_nb.groups()[0])
                        #if tr:
                            #print tr
                    if found_tr and found_nb:
                        if found_nb > genes[gene_id][1]:
                            genes[gene_id]=(found_tr, found_nb)
    gtf.close()
    transcripts = {}
    for item in genes:
        transcripts[genes[item][0]]=1
    return transcripts


def filter_gtf(trs_list, gtf_file, out_file):

    gtf = open(gtf_file, 'r')
    out = open(out_file, 'w')

    for line in gtf:

        split_line = re.split("\s+", line)
        if split_line[2] == "gene":
            out.write(line)
        elif line.startswith("#"):
            out.write(line)
        else:
            last_elt = ("").join(split_line[8:])
            attributes = re.split('\s*;\s*',last_elt)
            for attribute in attributes:
                #tr = re.search('^transcript_id"(\w+\.*\w*)"',attribute)
                tr = re.search('^transcript_id"(\w+\.*\w*)"',attribute)
                if tr:
                    if tr.groups()[0] in trs_list:
                        out.write(line)

    gtf.close()
    out.close()
    return




# =====================================================================================================================
#        MAIN
# =====================================================================================================================
parser = argparse.ArgumentParser()
parser.add_argument('-g', dest='gtf', required=True, help='gtf to filter', type=str)
parser.add_argument('-o', dest='output', required=True, help='output file', type=str)
options = parser.parse_args()

print("Selecting best transcript...")
transcripts = best_tr(options.gtf)
print("Filtering gtf...")
filter_gtf(transcripts,options.gtf,options.output)
