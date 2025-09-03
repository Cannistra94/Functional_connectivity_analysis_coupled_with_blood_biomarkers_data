#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

bet ${filename}.nii bet_${filename}.nii 

done
