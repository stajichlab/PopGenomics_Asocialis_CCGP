#script for visualizing admixture results
#Jan. 11th 2024
#jnadams

# Load libraries
library(tidyverse)
library(cowplot)

#   1. Read in the Q matrix data and tidy it such that each row is an observation 
#   (i.e., each row only has one individual and one K proportion). The easiest way 
#   to do so is with the pivot_longer() function within the tidyr package.
admix <- read.table("vcf/ASCCGP_v5.All.SNP.selected.4.Q", header = FALSE) %>%
  mutate(INDV = c(1:142))
admix_v2 <- read.table("vcf/ASCCGP_v5.All.SNP.selected.4.Q", header = FALSE) %>%
  mutate(INDV = admix_fam$V1, INDV = str_extract(INDV, "(\\d+$)"), INDV = case_when(is.na(INDV) ~ "55", TRUE ~ INDV))
admix_fam <- read.table("vcf/ASCCGP_v5.All.SNP.selected.fam", header = FALSE)
head(admix)
dim(admix)
head(admix_fam)
head(admix_v2)
tail(admix_v2)
View(admix_v2)

admix_tidy <- admix %>% 
  pivot_longer(names_to = "cluster", values_to = "proportion", -INDV)
head(admix_tidy)

admix_tidy_v2 <- admix_v2 %>% 
  pivot_longer(names_to = "cluster", values_to = "proportion", -INDV)
head(admix_tidy)
head(admix_tidy_v2)

#   2. Now, you should be able to build your plot using geom_bar with the 
#   position and stat args set to "stack" and "identity", respectively,
#   filling by cluster.
admix_tidy %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

admix_tidy_v2 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

admix_tidy %>%
  arrange(cluster, desc(proportion)) %>% 
  mutate(INDV = factor(INDV)) %>%
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

# For a cleaner/prettier structure plot:
admix_tidy %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.line=element_line(colour="black"),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # legend.position = "none",
        panel.border = element_rect(fill=NA, colour = "black", linetype = "solid", linewidth = 1.5),
        strip.text.y = element_text(size = 30, face = "bold"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.spacing = unit(-0.1, "lines"))

#code for making an admixture plot from a different tutorial
tbl=read.table("~/Data/MyProject/admixture/MyProject.HO.merged.6.Q")
indTable = read.table("~/Data/MyProject/admixture/MyProject.HO.merged.ind",
                      col.names = c("Sample", "Sex", "Pop"))
popGroups = read.table("~/Google Drive/GA_Workshop Jena/HO_popGroups.txt", col.names=c("Pop", "PopGroup"))

mergedAdmixtureTable = cbind(tbl, indTable)
mergedAdmWithPopGroups = merge(mergedAdmixtureTable, popGroups, by="Pop")
ordered = mergedAdmWithPopGroups[order(mergedAdmWithPopGroups$PopGroup),]
barplot(t(as.matrix(subset(ordered, select=V1:V6))), col=rainbow(6), border=NA)

#   1. Read in the Q matrix data and tidy it such that each row is an observation 
#   (i.e., each row only has one individual and one K proportion). The easiest way 
#   to do so is with the pivot_longer() function within the tidyr package.
admix_k2 <- read.table("vcf/ASCCGP_v5.All.SNP.selected.2.Q", header = FALSE) %>%
  mutate(INDV = c(1:142))
head(admix_k2)
dim(admix_k2)

admix_tidy_k2 <- admix_k2 %>% 
  pivot_longer(names_to = "cluster", values_to = "proportion", -INDV)
head(admix_tidy_k2)

#   2. Now, you should be able to build your plot using geom_bar with the 
#   position and stat args set to "stack" and "identity", respectively,
#   filling by cluster.
admix_tidy_k2 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

admix_tidy_k2 %>%
  arrange(cluster, desc(proportion)) %>% 
  mutate(INDV = factor(INDV)) %>%
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

# For a cleaner/prettier structure plot:
admix_tidy_k2 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.line=element_line(colour="black"),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # legend.position = "none",
        panel.border = element_rect(fill=NA, colour = "black", linetype = "solid", linewidth = 1.5),
        strip.text.y = element_text(size = 30, face = "bold"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.spacing = unit(-0.1, "lines"))

#   1. Read in the Q matrix data and tidy it such that each row is an observation 
#   (i.e., each row only has one individual and one K proportion). The easiest way 
#   to do so is with the pivot_longer() function within the tidyr package.
admix_k3 <- read.table("vcf/ASCCGP_v5.All.SNP.selected.3.Q", header = FALSE) %>%
  mutate(INDV = c(1:142))
head(admix_k3)
dim(admix_k3)

admix_tidy_k3 <- admix_k3 %>% 
  pivot_longer(names_to = "cluster", values_to = "proportion", -INDV)
head(admix_tidy_k3)

#   2. Now, you should be able to build your plot using geom_bar with the 
#   position and stat args set to "stack" and "identity", respectively,
#   filling by cluster.
admix_tidy_k3 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

admix_tidy_k3 %>%
  arrange(cluster, desc(proportion)) %>% 
  mutate(INDV = factor(INDV)) %>%
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

# For a cleaner/prettier structure plot:
admix_tidy_k3 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.line=element_line(colour="black"),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # legend.position = "none",
        panel.border = element_rect(fill=NA, colour = "black", linetype = "solid", linewidth = 1.5),
        strip.text.y = element_text(size = 30, face = "bold"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.spacing = unit(-0.1, "lines"))

#   1. Read in the Q matrix data and tidy it such that each row is an observation 
#   (i.e., each row only has one individual and one K proportion). The easiest way 
#   to do so is with the pivot_longer() function within the tidyr package.
admix_k5 <- read.table("vcf/ASCCGP_v5.All.SNP.selected.5.Q", header = FALSE) %>%
  mutate(INDV = c(1:142))
head(admix_k5)
dim(admix_k5)

admix_tidy_k5 <- admix_k5 %>% 
  pivot_longer(names_to = "cluster", values_to = "proportion", -INDV)
head(admix_tidy_k3)

#   2. Now, you should be able to build your plot using geom_bar with the 
#   position and stat args set to "stack" and "identity", respectively,
#   filling by cluster.
admix_tidy_k5 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

admix_tidy_k5 %>%
  arrange(cluster, desc(proportion)) %>% 
  mutate(INDV = factor(INDV)) %>%
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity")

# For a cleaner/prettier structure plot:
admix_tidy_k5 %>% 
  ggplot(aes(x = INDV, y = proportion, fill = cluster)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.line=element_line(colour="black"),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # legend.position = "none",
        panel.border = element_rect(fill=NA, colour = "black", linetype = "solid", linewidth = 1.5),
        strip.text.y = element_text(size = 30, face = "bold"),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.spacing = unit(-0.1, "lines"))
