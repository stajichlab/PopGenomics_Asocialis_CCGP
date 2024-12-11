#!/usr/bin/bash -l
#SBATCH -p intel -N 1 -n 16 --mem 32gb --out logs/make_vcf2phylip.plink.%a.log --time 120:00:00

MEM=32g


python vcf2phylip.py --input CCGP_output_filtered_pruned_v20.vcf --fasta --nexus --nexus-binary
