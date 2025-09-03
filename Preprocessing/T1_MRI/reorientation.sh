#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

fslreorient2std ${filename}/anat.nii reo_${filename}

done
