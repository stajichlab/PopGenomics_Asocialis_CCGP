#!/usr/bin/bash -l
#SBATCH -p short -c 2 -N 1 -n 1 --mem 4gb --out logs/index.log

module load samtools
module load bwa-mem2
if [ -f config.txt ]; then
	source config.txt
fi
grep ">" $REFGENOME | perl -p -e 's/^>(JAMYGW020+(\d+))\.\d+\s+/$1,$2\n/; s/>MT/MT,MT/' > $GENOMEFOLDER/chrom_nums.csv
echo "working off $REFGENOME - check if these don't match may need to update config/init script"

if [[ ! -f $REFGENOME.fai || $REFGENOME -nt $REFGENOME.fai ]]; then
	samtools faidx $REFGENOME
fi
if [[ ! -f $REFGENOME.bwt || $REFGENOME -nt $REFGENOME.bwt ]]; then
	bwa-mem2 index $REFGENOME
fi

DICT=$GENOMEFOLDER/$(basename $REFGENOME .fasta)".dict"
echo "$DICT"
if [[ ! -f $DICT || $REFGENOME -nt $DICT ]]; then
	rm -f $DICT
	samtools dict $REFGENOME > $DICT
	ln -s $(basename $DICT) $REFGENOME.dict
fi
