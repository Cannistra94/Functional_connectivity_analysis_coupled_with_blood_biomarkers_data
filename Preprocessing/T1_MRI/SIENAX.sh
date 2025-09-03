#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

sienax flirt_${filename}.nii.gz -r

done
