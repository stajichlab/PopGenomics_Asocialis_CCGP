#!/usr/bin/bash -l
#SBATCH -p intel -N 1 -n 16 --mem 32gb --out logs/make_vcf2phylip.plink_v2.%a.log --time 120:00:00

MEM=32g


python vcf2phylip.py --input 0003.vcf --fasta --nexus --nexus-binary
