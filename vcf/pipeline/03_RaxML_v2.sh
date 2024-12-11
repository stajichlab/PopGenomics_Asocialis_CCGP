#!/usr/bin/bash -l
#SBATCH -p stajichlab -N 1 -n 16 --mem 32gb --out logs/RaxML.%a.log --time 7-00:00:00

module load raxml/8.2.12
raxmlHPC-PTHREADS -T 10 -N 1000 -n myMLJob -s 142genomes_Jan20.phy -m MULTICAT -f a -x 12345 -p 12345
