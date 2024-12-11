#!/usr/bin/bash -l
#SBATCH -p short -c 2 -N 1 -n 1 --mem 4gb --out logs/admixture.log 

module load admixture/1.3.0

for K in $(seq 1 12); do admixture --cv ASCCGP_v5.All.SNP.selected.bed $K; done 
