#cross validation scores version 2
#Jan. 22nd 2024
#jnadams

#load require libs
library(stringr)
library(ggplot2)
library(dplyr)

#read in cross validation file
cv <- read.table("vcf/logs/cross_validation.txt")
head(cv)
tail(cv)

cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- select(cv, V4,K)
value = CV$V4
K <- str_sub(CV$K, 3)
K = as.numeric(K)
CV_new = tibble(value, K)
View(CV)
CV <- CV %>%
  mutate(K = str_sub(K, 3))
#K = as.numeric(K))
CV <- as.numeric(CV$K)
colnames(CV) <- c("CV","K")

graph_title="Cross-Validation plot"
x_title="K"
y_title="Cross-validation error"
graph_1<-ggplot(CV_new,aes(x=K,y=value))
graph_1
graph_1+geom_line()
graph_1+geom_line()+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))
