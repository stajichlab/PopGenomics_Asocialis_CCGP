#!/usr/bin/env Rscript

library(tidyverse)
library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(patchwork)
library(dplyr)
library(RColorBrewer)

getPalette1 = colorRampPalette(brewer.pal(12, "Paired"))

meta <- read_csv("CCGP_Herbarium_Metadata.csv", skip=1,col_names=c("ID","SampleID","Longitude","Latititude","Date","Locality",
                                                                  "General_Locality","NorthSouth","EastWest","Regions","Species"),
                 col_types = "icddDccccc")
# fix node IDs to match the tree names
CCGP = meta %>% filter(str_detect(SampleID, "^CCGP")) %>% mutate(label = gsub("CCGP", "JNA_AS_CCGP_",SampleID))
Herb = meta %>% filter(str_detect(SampleID, "^Herb")) %>% mutate(label = str_extract(SampleID, "[^_]+$"))
ref = meta %>% filter(str_detect(SampleID,"^ref") | str_detect(SampleID,"^Big") | str_detect(SampleID,"^Pink")) %>% mutate(label = SampleID)
meta = bind_rows(bind_rows(CCGP,Herb,ref)) %>% mutate(Sample_Locality = str_c(label," -> ",Locality))

rootname = "refJNA_AS_CCGP158"

#goodmappers
goodmappersTree <- read.newick("strain_tree/ASCCGP_v3.GoodMapper.SNP.fasttree.tre")

goodmappersTree <- root(goodmappersTree, outgroup = rootname, edgelabel = TRUE)

treeT <- left_join(as_tibble(goodmappersTree),meta,by="label")
treeT
treeT <- as.treedata(treeT)
treeT

#ggtree(treeT,color = "black", size = 1, linetype = 1) + geom_tiplab(fontface = "bold", size = 3, offset = 0.001) + xlim(0,1)

p1 <- ggtree(treeT) +
  geom_tippoint(
    mapping = aes(
      color = Regions,
      shape = EastWest,
    )) + scale_colour_manual(values=getPalette1(11)) +
      xlim(0, 1)  +
  scale_size_continuous(range = c(3, 8)) +
  geom_tiplab(fontface = "bold",  size = 2, offset = 0.001,
          mapping = aes(label = Sample_Locality)) +
  
  theme(legend.position = "right")  

#p1
ggsave('CCGP_tree_colors.pdf',p1,height=12,width=12)
#
#  
#
