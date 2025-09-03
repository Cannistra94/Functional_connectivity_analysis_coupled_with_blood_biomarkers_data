#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

flirt -in bet_${filename} -ref /usr/pubsw/packages/fsl/5.0.10/data/standard/MNI152_T1_2mm_brain -out flirt_${filename} -omat flirt_${filename}.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

done
