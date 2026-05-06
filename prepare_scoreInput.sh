#!/bin/sh


jobid=$1
TF=$2



cd Output_$jobid/threads_$jobid


mkdir $TF"_threads"
mkdir $TF"_scores"
cd $TF"_threads"

threadfile=$(ls ../*${TF}*txt | head -1)
pdbfile=$(ls ../*${TF}*pdb | head -1)
cp $threadfile .
cp $pdbfile .

python ../../../prepare_scoreInput.py $TF $threadfile $jobid

