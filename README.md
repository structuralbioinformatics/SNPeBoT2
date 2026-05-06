#     SNPeBoT2: SNPeBoT2 to be deployed as a downloadable Package      

## Installation instructions:
1) download this repository
2) install ModCRElib within the SNPeBoT2 folder (follow ModCRElib installation instructions from https://github.com/structuralbioinformatics/ModCRElib)
	- In the current directory you will have a folder labeled ModCRElib that contains everything from the github) 
	- Ensure all ModCRElib dependencies are met
3) in the AssociationFilter folder download and decompress the following:
	- https://cdn.netbiol.org/tflink/download_files/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz
	- https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz

## Running SNPeBoT2
Predict_SNP_effect.sh:
	This is the bash script that will generate ASB predictions on a list of reference and alternate sequence for a given Transcription Factor
	It takes an input of: 
		1) TF name: This is the name or id of the transcription factor, it should match the id or name used in 3 and 4
		2) job id: A unique job id that will be used to designate output folder names 
		3) ASB file location: The path to the file containing the SNPs to be tested.
			- file has 4 columns per SNP: 1) TF name 
			                              2) reference sequence (51 nt long with mutation position at 26) 
			                              3) alternate sequence (51 nt long with mutation position at 26) 
			                              4) SNP id (suggest {chr}-{position}-{ref allele}-{alt allele}-{TF name}
FIGLA   GGATTTATCCATGTTTTTGCATGTAGAGATGGCTTGAAAAACAAACTATAT     GGATTTATCCATGTTTTTGCATGTACAGATGGCTTGAAAAACAAACTATAT     chr7-133924858-G-C-FIGLA 
		4) model of TF bound to DNA: This can be generated either with ModCRElib or another tool 
		5) id of template used in model: If ModCRElib was used this can be found in the model file name otherwise this must be designated "general".
	The models which must be submitted by the user will be stored automatically in Output_{job id}/models_{job id}
	The prediction of the effect on the binding of the TF that each SNP will have are stored in Output_{job id}/Predictions_{job id}.tsv 
	

RunAssoc.sh:
	This is the script that will predict if the SNP is associated with a disease.
	It takes an input of:
		1) TF name: Use uniprot id of the transcription factor
		2) Chrom: the chromosome the SNP is located on 
		3) Prediction: the prediction given by Predict_SNP_effect.sh (gain,loss,no-change)
		4) Position: The position of the SNP
		5) jobid: The job id submitted in Predict_SNP_effect.sh
		6) assembly: The human genome assembly used for the SNP location
		7) SNPId: The id for the SNP (as used in Predict_SNP_effect.sh)
	The output will be stored in Output_$jobid/associations/assoc_$SNPId.csv




	disgenet/all_gene_disease_pmid_associations.tsv.gz


regulator info from:
	https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000025/QTD000434/QTD000434.credible_sets.tsv.gz

