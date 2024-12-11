#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 16 --mem 32gb --out logs/bcftools%a.log --time 168:00:00

module load bcftools/1.18

bcftools isec CCGP_output_filtered_pruned_v20.vcf.gz Astr_calls.vcf.gz -p Astr_overlap # overlap Astr to Asoc
bcftools view -Oz -o Astr_in_Asoc.vcf.gz Astr_overlap/0003.vcf # save the set in common 

bcftools index Astr_in_Asoc.vcf.gz # index
bcftools isec Astr_in_Asoc.vcf.gz CCGP_output_filtered_pruned_v20.vcf.gz -p Astr_overlap_backwards # get the Asoc strain subset of sites also in Astr
