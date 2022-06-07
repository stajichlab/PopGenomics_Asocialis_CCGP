#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 96 --mem 384gb --out logs/mmseqs_classify_reads.%a.log -a 1-28

module load workspace/scratch
module load mmseqs2
module load KronaTools

if [ -f config.txt ]; then
  source config.txt
fi

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi
INDIR=input
DB2=/srv/projects/db/ncbi/mmseqs/uniref50
DB2NAME=$(basename $DB2)
IFS=,
OUTSEARCH=results/mmseqs2
mkdir -p $OUTSEARCH
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
    mkdir -p $OUTSEARCH/$STRAIN
    BASEPATTERN=$(echo $FILEBASE | perl -p -e 's/\;/ /g; ')
    if [ ! -s $OUTSEARCH/$STRAIN/mmseq_${DB2NAME}_report ]; then
	mmseqs easy-taxonomy $INDIR/$BASEPATTERN $DB2 $OUTSEARCH/$STRAIN/mmseq_$DB2NAME $SCRATCH --threads $CPU --lca-ranks kingdom,phylum,family  --tax-lineage 1
	ktImportTaxonomy -o $OUTSEARCH/$STRAIN/mmseq_${DB2NAME}.krona.html $OUTSEARCH/$STRAIN/mmseq_${DB2NAME}_report
    fi
done
