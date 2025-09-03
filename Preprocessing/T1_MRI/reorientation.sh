#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

fslreorient2std ../../T1/Knee_Pain/${filename}/anat.nii /autofs/space/yintang_003/users/valeria/Lab_cohorts/derivatives/T1/reo_testoldversion/reo_${filename}

done
