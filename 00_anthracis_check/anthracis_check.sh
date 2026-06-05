#!/bin/bash

path='/mnt/DATA/eurothrax/01_nxtthrax/results/02_anthracis_check'

mkdir -p ${path}/blast_results

for genome in $(ls ${path}/genomes/*.fna); do
	
	bn=$(basename ${path}"$genome" .fna)
	
	echo "Working on:" ${bn}
	
	echo ${bn} >> dhp61_results.txt
	
	blastn -db ./dhp61/dhp61 -query $genome -out ${path}/blast_results/${bn}_dhp61.out
	
	grep Identities ${path}/blast_results/${bn}_dhp61.out >> dhp61_results.txt
	
	echo ${bn} >> pl3_results.txt
	
	blastn -db ./pl3/pl3 -query $genome -out ${path}/blast_results/${bn}_pl3.out
		
	grep Identities ${path}/blast_results/${bn}_pl3.out >> pl3_results.txt
done


