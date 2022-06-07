#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 8 --mem 64gb --out logs/bracken_unmapped.%a.log -a 1-159

module load workspace/scratch
module load bracken
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
DB=pluspf
INDIR=unmapped
IFS=,
OUTDIR=results/bracken_unmapped
READ_LEN=150
#CLASSIFICATION_LEVEL=G
THRESHOLD=10
mkdir -p $OUTDIR
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
    if [ ! -s $OUTDIR/$STRAIN.kreport2 ]; then
	kraken2 -d $DB --threads $CPU --report  $OUTDIR/$STRAIN.kreport2 $INDIR/$STRAIN.fastq.gz $INDIR/${STRAIN}_single.fastq.gz > $OUTDIR/$STRAIN.kraken2 
    fi
    ktImportTaxonomy -o $OUTDIR/$STRAIN.kraken2.krona.html $OUTDIR/$STRAIN.kraken2

    for CLASSIFICATION_LEVEL in G S;
    do
	if [ ! -f $OUTDIR/${STRAIN}.${CLASSIFICATION_LEVEL}.bracken ]; then
	    time bracken -d $KRAKEN2_DB_PATH/${DB} -i $OUTDIR/${STRAIN}.kreport2 -o $OUTDIR/${STRAIN}.${CLASSIFICATION_LEVEL}.bracken -r ${READ_LEN} -l ${CLASSIFICATION_LEVEL} -t ${THRESHOLD}
	fi
    done
done
