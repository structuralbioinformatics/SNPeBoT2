#!bin/bash

TF=$1
Chrom=$2
Prediction=$3
Position=$4
jobid=$5
assembly=$6
SNPId=$7
echo "python AssociationFilter/AssociationFilter.py $TF $Chrom $Prediction $Position $assembly Output_$jobid/associations/assoc_$SNPId.csv"

mkdir Output_$jobid/associations
python AssociationFilter/AssociationFilter.py $TF $Chrom $Prediction $Position $assembly Output_$jobid/associations/assoc_$SNPId.csv

