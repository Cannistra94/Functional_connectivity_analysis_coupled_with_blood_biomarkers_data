#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

echo ${filename},`fslstats ${filename}/${filename}_pve_2 -M -V | awk '{print $1 * $3}'` >> WM.csv

done
