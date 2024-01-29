# clean environment variables
rm(list=ls()) 
library(readxl)
library(ggplot2)
library("dplyr")

data1 <- read.table("The differential KEGG pathways_LDA.txt",sep="\t",header = T)

data_Metabolism<- data1[data1$L1=="Metabolism",]


colnames(data_Metabolism)[3]<-"SamplingSeason"
data_Metabolism$SamplingSeason <- factor(data_Metabolism$SamplingSeason,levels = c("Spring", "Autumn")) 

data_Metabolism<-data_Metabolism %>% 
  arrange(data_Metabolism$`The value of log10 with the greatest average abundance`)

level <- as.matrix(data_Metabolism[,1])  
level <- as.vector(level)

data_Metabolism$KEGG <- factor(data_Metabolism$KEGG,levels =level)

p<- ggplot(data_Metabolism,aes(x=SamplingSeason,y=KEGG,size=abundance,color=L2))  + 
  geom_point(shape=16) + theme_bw()+
  scale_size(range = c(4, 10))
  
p<-p
p

#ggsave("Metabolism_baboo.pdf", p, width = 10, height = 6)

