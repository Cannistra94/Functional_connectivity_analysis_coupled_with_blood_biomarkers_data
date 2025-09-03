#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

run_first_all -i flirt_${filename}.nii -o ${filename}.nii -b

done
