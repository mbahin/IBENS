#!/bin/bash

tot = 0

echo "Unknown / Unknown"
n=$(sed 1d $1 | awk '$2 == "Unknown" && $3 == "Unknown"' | wc -l)
echo $n
tot=$(($tot + $n))
echo "Unknown / UMI"
n=$(sed 1d $1 | awk '$2 == "Unknown" && $3 != "Unknown"' | wc -l)
echo $n
tot=$(($tot + $n))
echo "Trimmed / Unknown"
n=$(sed 1d $1 | awk '$2 == "trimmed" && $3 == "Unknown"' | wc -l)
echo $n
tot=$(($tot + $n))
echo "Trimmed / UMI"
n=$(sed 1d $1 | awk '$2 == "trimmed" && $3 != "Unknown"' | wc -l)
echo $n
tot=$(($tot + $n))
echo "UMI / Unknown"
n=$(sed 1d $1 | awk '$2 != "trimmed" && $2 != "Unknown" && $3 == "Unknown"' | wc -l)
echo $n
tot=$(($tot + $n))
echo "UMI / UMI (same)"
n=$(sed 1d $1 | awk '$2 != "trimmed" && $2 != "Unknown" && $3 != "Unknown" && $3 == $2' | wc -l)
echo $n
tot=$(($tot + $n))
echo "UMI / UMI (diff)"
n=$(sed 1d $1 | awk '$2 != "trimmed" && $2 != "Unknown" && $3 != "Unknown" && $3 != $2' | wc -l)
echo $n
tot=$(($tot + $n))
echo "Total: $tot"
