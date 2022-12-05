#!/usr/bin/bash -l
#sbatch -p short -N 1 -n 24 --mem 96gb --out logs/cram_stats.log

CPU=24
module load samtools
module load parallel

parallel -j $CPU samtools flagstat {} \> {}.stats.txt ::: $(ls aln/*.cram)
