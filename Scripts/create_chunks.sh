#!/bin/bash

# Arguments
# $1 -> input file to cut
# $2 -> chunk size (line count)
# $3 -> output file prefix
# $4 -> output file suffix

# Getting the total file lines
nb_lines=$(wc -l $1 | awk '{print $1}')

# Initializing the file index to 0
n=0

# Creating the file chunks
for i in $(seq $2 $2 $nb_lines); do
	head -n $i $1 | tail -n $2 > $3.$n.$4
	n=$(($n + 1))
done

# Creating the last chunk if necessary
if [[ $(($nb_lines % $2)) != 0 ]]; then
	tail -n $(($nb_lines % $2)) $1 > $3.$n.$4
fi
