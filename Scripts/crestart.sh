#!/bin/bash

echo "### General statistics ###"
condor_q -g -allusers -l | grep -i numjobstart | sort | uniq -c

echo -e "\n### Details on restarted jobs ###"
cqga -l | egrep -i "^owner|^numjobstart" | awk '{if(NR%2){line=$0}else{print line"\t"$0}}' | sed 's/NumJobStarts = \([0-9]*\)\tOwner = "\([^"]*\)"/\2\t\1/g' | awk '$2 > 1' | sort | uniq -c | awk '{print $2"\t"$3 " times\tcount="$1}' | column -t
