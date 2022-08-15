#!/bin/bash

#SBATCH --partition=batch
#SBATCH --cpus-per-task=24
#SBATCH --mem=200g
#SBATCH --time=0-72:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="kracken"

# Building custom Kraken database
module load kraken2/2.1.2

DBNAME="custom_DB"
FNA="/bigdata/stajichlab/shared/projects/Lichen/Acarospora_socialis_HerbariumSamples/custom_DB/library/lichen/lichen.fna"

kraken2-build --add-to-library $FNA --db $DBNAME

kraken2-build --build --threads 24 --db $DBNAME

echo "custom build complete"
