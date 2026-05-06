#!/usr/bin/env bash
set -euo pipefail

TF=$1
jobid=$2
SequenceFile=$3
pdb=$4
template=$5

echo "Running ModCRElib pipeline for $TF on $jobid"
mkdir "Output_"$jobid
mkdir "Output_"$jobid"/models_"$jobid
mv $pdb "Output_"$jobid"/models_"$jobid"/"
pdb="Output_${jobid}/models_${jobid}/${pdb}"


#TF=AHR
#jobid=1866
#SequenceFile="Resources/AHR_ASBs.in"
echo "Generating PWMs from Model"
cd ModCRElib
echo "python RenamePDB.py $pdb $template $TF"
newdir=$(python RenamePDB.py $pdb $template $TF)
echo $newdir $pdb
rm "../$pdb"
pdb=$newdir
if [[ "$template" == "general" ]]; then
    echo "[INFO] Template is 'general' — running pwmBoilerGeneral.sh"
    bash pwmBoilerGeneral.sh "$jobid"
    echo "PWMs located in Output_$jobid/pwm_$jobid"
    echo "Generating thread files"
    bash get_best_bindings_threadsBoilerGeneral.sh $jobid
    echo "thread file located in Output_$jobid/threads_$jobid"
else
    echo "[INFO] Template is '$template' — running pwmBoiler.sh"
    bash pwmBoiler.sh "$jobid"
    echo "PWMs located in Output_$jobid/pwm_$jobid"
    echo "Generating thread files"
    bash get_best_bindings_threadsBoiler.sh $jobid
    echo "thread file located in Output_$jobid/threads_$jobid"
fi
#bash pwmBoiler.sh $jobid
#echo "PWMs located in Output_$jobid/pwm_$jobid"
#echo "Generating thread files"
#bash get_best_bindings_threadsBoiler.sh $jobid
#echo "thread file located in Output_$jobid/threads_$jobid"
cd ..
#END of Model Handeling

#Now getting all of necessary binding site sequences


threadfiles=( $(ls Output_$jobid/threads_$jobid/*txt) )
threadfile="${threadfiles[0]}"
BindingSequence="${threadfile##*_}"
BindingSequence="${BindingSequence%.txt}"
BindingLength=${#BindingSequence}
echo "Finished first stage"
python Getinput.py --ASBInput $SequenceFile --tf $TF --length $BindingLength --jobID $jobid



bash prepare_scoreInput.sh ${jobid} "$TF"
echo "threads : Output_$jobid/threads_$jobid/${TF}_threads      scores : Output_$jobid/threads_$jobid/${TF}_scores" 

rm -f Output_$jobid/threads_$jobid/${TF}_threads/*txt

#exit 0

cd Output_$jobid/threads_$jobid/${TF}_threads || exit 1
python create_threads.py
total=$(($(ls | wc -l)-2))
threadfiles=( $(ls *txt) )
cd ../../..
echo "Finished prep for $TF"
echo "total thread files created: $total"


cp ScorerBoiler.sh  Scorer_${jobid}.sh

echo "dummy=\"${jobid}_dummy\"" >> "Scorer_${jobid}.sh"

if [[ "$template" == "general" ]]; then
	for thread in "${threadfiles[@]}"; do
		echo "python ModCRElib/exe/scorer.py --threading --dummy \$dummy -i Output_$jobid/threads_$jobid/${TF}_threads/$thread  --norm -o Output_$jobid/threads_$jobid/${TF}_scores/${thread%.txt} --pdb \$pdb --pbm \$pbm  --verbose " >> "Scorer_${jobid}.sh"
	done
else
	for thread in "${threadfiles[@]}"; do
                echo "python ModCRElib/exe/scorer.py --threading --dummy \$dummy -i Output_$jobid/threads_$jobid/${TF}_threads/$thread  --norm -o Output_$jobid/threads_$jobid/${TF}_scores/${thread%.txt} --pdb \$pdb --pbm \$pbm  --verbose  --auto --known " >> "Scorer_${jobid}.sh"
	done
fi

#for thread in "${threadfiles[@]}"; do
#    echo "python ModCRElib/exe/scorer.py --threading --dummy \$dummy -i Output_$jobid/threads_$jobid/${TF}_threads/$thread  --norm -o Output_$jobid/threads_$jobid/${TF}_scores/${thread%.txt} --pdb \$pdb --pbm \$pbm  --verbose  --auto --known " >> "Scorer_${jobid}.sh"
#done



bash Scorer_${jobid}.sh

firstfile="${threadfiles[0]}"
prefix="${firstfile%_*}_" 

python processScores.py $TF $prefix $jobid
python predict.py  Output_$jobid/ModelInput_$jobid.tsv  Output_$jobid/Predictions_$jobid.tsv

rm -rf ${jobid}_dummy_[0-9]*
rm -rf Scorer_${jobid}.sh


