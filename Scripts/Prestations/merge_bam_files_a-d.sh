#! /bin/bash

# Mathieu Bahin, 08/09/16
# Script to merge the a, b, c and d files from the sample produced by the genomics platform.

# Getting commandline options back
while getopts "d:" OPTION
do
        case $OPTION in
        d) input_dir=$OPTARG;;
    esac
done

# Merging the BAM files for each sample
for sample in $(ls $input_dir/*bam | sed 's/.*\/filtered_mapper_results_\(.*\)\-[abcd].bam/\1/g' | sort | uniq); do
	echo 'Processing sample '$sample'...'
	samtools merge -cp $sample.merged.bam $input_dir/*$sample*bam
done
