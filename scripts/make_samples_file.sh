#!/usr/bin/bash -l
echo "STRAIN,FILEBASE" > samples.csv
pushd input
#ls *R1.fastq.gz | perl -p -e 's/((\S+)_R1\.fastq\.gz)/$2,$2_R[12].fastq.gz/' >> ../samples.csv
ls *R1_001.fastq.gz | perl -p -e 's/(((\S+)_S\d+)_L\d+_R1_001\.fastq\.gz)/$3,$2\_R[12]_001.fastq.gz/'  >> ../samples.csv
popd
perl scripts/fix_CCGP_multirun.pl samples.csv > sample_fix.csv
