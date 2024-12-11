#!/usr/bin/bash -l
#SBATCH -p epyc -c 8 -N 1 -n 1 --mem 96gb --out logs/algatr_run.%A.log

Rscript Algatr_test_v2.R
