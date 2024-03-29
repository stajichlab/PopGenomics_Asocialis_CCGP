#!/usr/bin/bash -l
#SBATCH --nodes 1 --ntasks 24 --time 2:00:00 -p short --mem 64G --out logs/mosdepth.parallel.log
#SBATCH -J mosdepth

# This one goes to 11
module load parallel

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
 CPU=2
fi
module load mosdepth
mkdir -p coverage/mosdepth
source config.txt
if [ 0 ]; then
for WINDOW in 5000 10000 50000
do
	parallel --jobs $CPU mosdepth -f $REFGENOME -T 1,10,50,100,200 -n --by $WINDOW -t 2 "{= s:$ALNFOLDER\/:coverage/mosdepth/:; s:\.$HTCEXT:.${WINDOW}bp: =}" {} ::: $ALNFOLDER/*.$HTCEXT
done
fi
mkdir -p plots
Rscript scripts/plot_mosdepth_CNV.R

