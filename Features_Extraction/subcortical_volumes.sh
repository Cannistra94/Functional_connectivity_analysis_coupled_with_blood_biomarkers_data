#!/bin/bash

# read inputs

list=$1

dos2unix $list

for filename in `cat $list`; do

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 25.5 -u 26.5 -V | awk '{print $2}'` >> Vol_LAccu.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 57.5 -u 58.5 -V | awk '{print $2}'` >> Vol_RAccu.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 17.5 -u 18.5 -V | awk '{print $2}'` >> Vol_LAmyg.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 53.5 -u 54.5 -V | awk '{print $2}'` >> Vol_RAmyg.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 15.5 -u 16.5 -V | awk '{print $2}'` >> Vol_Brain-Stem.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 10.5 -u 11.5 -V | awk '{print $2}'` >> Vol_LCaudate.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 49.5 -u 50.5 -V | awk '{print $2}'` >> Vol_RCaudate.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 16.5 -u 17.5 -V | awk '{print $2}'` >> Vol_LHippo.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 52.5 -u 53.5 -V | awk '{print $2}'` >> Vol_RHippo.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 12.5 -u 13.5 -V | awk '{print $2}'` >> Vol_LPallidum.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 51.5 -u 52.5 -V | awk '{print $2}'` >> Vol_RPallidum.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 12.5 -u 13.5 -V | awk '{print $2}'` >> Vol_LPutamen.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 50.5 -u 51.5 -V | awk '{print $2}'` >> Vol_RPutamen.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 9.5 -u 10.5 -V | awk '{print $2}'` >> LThal.csv

echo ${filename},`fslstats ${filename}_all_fast_firstseg -l 48.5 -u 49.5 -V | awk '{print $2}'` >> RThal.csv

done
