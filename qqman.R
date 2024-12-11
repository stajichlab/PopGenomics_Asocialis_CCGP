#input for qqman and qqman code
#Jan. 29th 2024
#jnadams

#get current working directory
getwd()
#read in fst file
socal_coastal_vs_desert <- read_tsv("vcf_popgenome/coastal_socal_vs_desert_socal.weir.fst", col_names = TRUE)
head(socal_coastal_vs_desert)

#replace NaN with 0
socal_coastal_vs_desert$WEIR_AND_COCKERHAM_FST[is.nan(socal_coastal_vs_desert$WEIR_AND_COCKERHAM_FST)]<-0
is.nan(socal_coastal_vs_desert$WEIR_AND_COCKERHAM_FST)
head(socal_coastal_vs_desert)
tail(socal_coastal_vs_desert)

