#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 8 --mem 32gb --out logs/unmapped_asm_viralsearch.%a.log -a 1-187

module load diamond
module load minimap2

MEM=32
UNMAPPEDASM=unmapped_asm
OUTSEARCH=unmapped_asm_search
mkdir -p $OUTSEARCH
DB=/bigdata/stajichlab/shared/lib/funannotate_db/uniprot.dmnd
VIRALDB=/bigdata/stajichlab/shared/lib/Viral/RefSeq/2023_05_05/viral.protein.dmnd
VIRALDNA=/bigdata/stajichlab/shared/lib/Viral/RefSeq/2023_05_05/viral.mmi
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

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE 
do
  if [ ! -f $OUTSEARCH/$STRAIN.uniprot.diamond.tsv ]; then
    diamond blastx --db $DB -q $UNMAPPEDASM/$STRAIN/scaffolds.fasta --out $OUTSEARCH/$STRAIN.uniprot.diamond.tsv \
  -f 6 -b12 -c1 --memory-limit $MEM --ultra-sensitive --long-reads --threads $CPU
  fi
  if [ ! -f $OUTSEARCH/$STRAIN.RefGenome.paf ]; then
    minimap2 $VIRALDNA $UNMAPPEDASM/$STRAIN/scaffolds.fasta -t $CPU --cs=long > $OUTSEARCH/$STRAIN.RefViral.minimap.paf
  fi
  if [ ! -f $OUTSEARCH/$STRAIN.Viral.diamond.tsv ]; then
	  diamond blastx --db $VIRALDB -q $UNMAPPEDASM/$STRAIN/scaffolds.fasta --out $OUTSEARCH/$STRAIN.Viral.diamond.tsv -f 6 -b12 -c1 --memory-limit $MEM --ultra-sensitive --long-reads --threads $CPU
  fi
done
