
# The relationship between microbial and diet beta diversity. 

# clean environment variables
rm(list=ls()) 

#16s
# load related packages
if(!require(ggplot2))install.packages("ggplot2")
library("ggplot2") 
library("dplyr")
library("ggpmisc")

#  Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#bray_curtis distance
betadiv <-  read.table("bray_curtis.txt", header=T, row.names= 1, sep="\t")


betadiv = betadiv[rownames(design), rownames(design)]

betadis <- as.matrix(betadiv)

betadis[upper.tri(betadis,diag = T)] <- NA
rownames(betadis)[rownames(betadis) == "A195_2"] <- "A1952"
colnames(betadis)[colnames(betadis) == "A195_2"] <- "A1952"


library(reshape2)

p <- melt(betadis)

frame <- na.omit(p)
#coi################################################################################################################################
 
library(vegan) 
library(tidyr)


otu_Flattening<-read.table("otu_rare_coi.txt", header=T, row.names= 1, sep="\t") 
 otu_Flattening<-t(otu_Flattening)


betadiv_coi <- vegdist(otu_Flattening,method = "bray")
betadiv_coi<- as.matrix(betadiv_coi)
design_coi <-  design[row.names(design)%in%row.names(betadiv_coi),] 



betadiv_coi = betadiv_coi[rownames(design_coi), rownames(design_coi)]
betadis_coi <- as.matrix(betadiv_coi)


betadis_coi[upper.tri(betadis_coi,diag = T)] <- NA
rownames(betadis_coi)[rownames(betadis_coi) == "A195_2"] <- "A1952"
colnames(betadis_coi)[colnames(betadis_coi) == "A195_2"] <- "A1952"


p <- melt(betadis_coi)
frame_coi <- na.omit(p)

colnames(frame_coi) <- c("Var1","Var2","value_coi")


frame_coi<-unite(frame_coi,"column1",c("Var1","Var2"), sep="-", remove = T)
frame<-unite(frame,"column1",c("Var1","Var2"), sep="-", remove = T)
frame_fan<-merge(frame_coi,frame,by="column1")

frame_fan<-frame_fan %>%
  mutate(SamplingSeason = case_when(
    grepl("^A\\d+-B\\d+$", column1) ~ "twoSeason",
    grepl("^A\\d+-A\\d+$", column1) ~ "Spring",
    grepl("^B\\d+-B\\d+$",column1) ~ "Autumn",
    grepl("^B\\d+-A\\d+$", column1) ~ "twoSeason",
    
  ))

frame_fan$SamplingSeason <- factor(frame_fan$SamplingSeason,levels = c("Spring", "Autumn","twoSeason"))



betadis_ma<- betadiv[row.names(betadiv)%in%row.names(betadiv_coi),colnames(betadiv)%in%colnames(betadiv_coi)]



frame_fan$allsample<- "allsample"
library("dplyr")
frame_fan<-frame_fan %>%
  mutate(group = case_when(
    grepl("Spring", SamplingSeason) ~ "Apre-breeding",
    grepl("Autumn",SamplingSeason) ~ "post-breeding",))
frame_fan$group <- factor(frame_fan$group,levels = c("Apre-breeding", "post-breeding"))  


p4<- ggplot(frame_fan,aes(x = value_coi, y = value))+
  
  geom_point(aes(color=group))+
  
  geom_smooth(frame_fan[frame_fan$SamplingSeason=="Spring",], mapping=aes(x = value_coi, y = value,color=group),method = "lm")+
  geom_smooth(frame_fan[frame_fan$SamplingSeason=="Autumn",], mapping=aes(x = value_coi, y = value,color=group),method = "lm")+
  geom_smooth(frame_fan, mapping=aes(x = value_coi, y = value,color=allsample),method = "lm")+
  scale_y_continuous(name = "Bray-Curtis microbiome dissimilarity")+
  scale_x_continuous(name = "Bray-Curtis dietary(coi) dissimilarity")
p4<-p4+
  scale_color_manual(values=c("blue","#FFCC00","#33CC33"))+
  annotate("text",label="pearson:r==0.1512~p<0.05",family="Times",parse=T,x=0.70,y=0.2,color="blue",size=4)
  
p4


#animal diet vs microbiota 
set.seed(000)
mantel(betadis_ma,betadiv_coi,method="pearson")#0.1521
set.seed(000)
mantel(betadis_ma[c(1:28),c(1:28)],betadiv_coi[c(1:28),c(1:28)],method="pearson")
set.seed(000)
mantel(betadis_ma[c(29:56),c(29:56)],betadiv_coi[c(29:56),c(29:56)],method="pearson")
######################################################rbcl######################################################
 
#bray_curtis distance
otu_Flattening<-read.table("otu_rare_rbcl.txt", header=T, row.names= 1, sep="\t") 
otu_Flattening<-t(otu_Flattening)

betadiv_rbcl <- vegdist(otu_Flattening,method = "bray")
betadiv_rbcl<- as.matrix(betadiv_rbcl)
design_rbcl <-  design[row.names(design)%in%row.names(betadiv_rbcl),] 



betadiv_rbcl = betadiv_rbcl[rownames(design_rbcl), rownames(design_rbcl)]
betadis_rbcl <- as.matrix(betadiv_rbcl)



betadiv_rbcl = betadiv_rbcl[rownames(design_rbcl), rownames(design_rbcl)]
betadis_rbcl <- as.matrix(betadiv_rbcl)

betadis_rbcl[upper.tri(betadis_rbcl,diag = T)] <- NA

rownames(betadis_rbcl)[rownames(betadis_rbcl) == "A195_2"] <- "A1952"
colnames(betadis_rbcl)[colnames(betadis_rbcl) == "A195_2"] <- "A1952"


p_rbcl <- melt(betadis_rbcl)
frame_rbcl <- na.omit(p_rbcl)

colnames(frame_rbcl) <- c("Var1","Var2","value_rbcl")


frame_rbcl<-unite(frame_rbcl,"column1",c("Var1","Var2"), sep="-", remove = T)

frame_fan_rbcl<-merge(frame_rbcl,frame,by="column1")

frame_fan_rbcl<-frame_fan_rbcl %>%
  mutate(SamplingSeason = case_when(
    grepl("^A\\d+-B\\d+$", column1) ~ "twoSeason",
    grepl("^A\\d+-A\\d+$", column1) ~ "Spring",
    grepl("^B\\d+-B\\d+$",column1) ~ "Autumn",
    grepl("^B\\d+-A\\d+$", column1) ~ "twoSeason",
    
  ))

frame_fan_rbcl$SamplingSeason <- factor(frame_fan_rbcl$SamplingSeason,levels = c("Spring", "Autumn","twoSeason"))


betadis_ma_rbcl<- betadiv[row.names(betadiv)%in%row.names(betadiv_rbcl),colnames(betadiv)%in%colnames(betadiv_rbcl)]


betadiv_coi<-as.data.frame(betadiv_coi)
betadiv_rbcl<-as.data.frame(betadiv_rbcl)



frame_fan_rbcl$allsample<-"allsample"

library("dplyr")
frame_fan_rbcl<-frame_fan_rbcl %>%
  mutate(group = case_when(
    grepl("Spring", SamplingSeason) ~ "Apre-breeding",
    grepl("Autumn",SamplingSeason) ~ "post-breeding",))




frame_fan_rbcl$group <- factor(frame_fan_rbcl$group,levels = c("Apre-breeding", "post-breeding"))  


p8<- ggplot(frame_fan_rbcl,aes(x = value_rbcl, y = value))+
  
  geom_point(aes(color=group))+
  
  geom_smooth(frame_fan_rbcl[frame_fan_rbcl$SamplingSeason=="Spring",], mapping=aes(x = value_rbcl, y = value,color=group),method = "lm")+
  geom_smooth(frame_fan_rbcl[frame_fan_rbcl$SamplingSeason=="Autumn",], mapping=aes(x = value_rbcl, y = value,color=group),method = "lm")+
  geom_smooth(frame_fan_rbcl, mapping=aes(x = value_rbcl, y = value,color=allsample),method = "lm")+
  scale_y_continuous(name = "Bray-Curtis microbiome dissimilarity")+
  scale_x_continuous(name = "Bray-Curtis dietary(rbcl) dissimilarity")
  

p8<-p8+
  scale_color_manual(values=c("blue","#FFCC00","#33CC33","#CCCCCC"))+
  annotate("text",label="pearson:r==0.2924~p==0.001",family="Times",parse=T,x=0.70,y=0.2,color="blue",size=4)



p8



library(patchwork)


p4<- p4+ ggtitle('A')
p8<- p8+ ggtitle('B')
p9<-p4+p8
p9
#ggsave("pearson_BC_diet_micro2024.1.18.pdf", p9, width =10, height = 4)
#save.image(file = "Figure4A,B.RData")
#plant vs microbial
set.seed(000)
mantel(betadis_ma_rbcl,betadiv_rbcl,method="pearson")
set.seed(000)
mantel(betadis_ma_rbcl[c(1:28),c(1:28)],betadiv_rbcl[c(1:28),c(1:28)],method="pearson")
set.seed(000)
mantel(betadis_ma_rbcl[c(29:58),c(29:58)],betadiv_rbcl[c(29:58),c(29:58)],method="pearson")

#animal vs plant
betadiv_rbcl<- betadiv_rbcl[row.names(betadiv_rbcl)%in%row.names(betadiv_coi),colnames(betadiv_rbcl)%in%colnames(betadiv_coi)]
set.seed(000)
mantel(betadiv_coi,betadiv_rbcl,method="pearson")
set.seed(000)
mantel(betadiv_coi[c(1:28),c(1:28)],betadiv_rbcl[c(1:28),c(1:28)],method="pearson")
set.seed(000)
mantel(betadiv_coi[c(29:56),c(29:56)],betadiv_rbcl[c(29:56),c(29:56)],method="pearson")


#  partial correlation  test
head(frame)
frame_animal_plant<-merge(frame_coi,frame_rbcl,by="column1")

frame_all<-merge(frame,frame_animal_plant,by="column1")

if(!require(ppcor))install.packages("ppcor")
library(ppcor)

colnames(frame_all)<-c("Samply","Microbiome","animal_diet","Plant_diet")


data2<-pcor(frame_all[,c(2,3,4)])
data<-data2$estimate
data.p<-data2$p.value

#season

frame_all<-frame_all %>%
  mutate(SamplingSeason = case_when(
    grepl("^A\\d+-B\\d+$", Samply) ~ "twoSeason",
    grepl("^A\\d+-A\\d+$", Samply) ~ "pre-breeding",
    grepl("^B\\d+-B\\d+$",Samply) ~ "post-breeding",
    grepl("^B\\d+-A\\d+$", Samply) ~ "twoSeason",))


#pre-breeding
pca_all_spring<-frame_all[frame_all$SamplingSeason=="pre-breeding",]
data2<-pcor(pca_all_spring[,c(2,3,4)])
data<-data2$estimate
data.p<-data2$p.value

#post-breeding
pca_all_autumn<-frame_all[frame_all$SamplingSeason=="post-breeding",]
data2<-pcor(pca_all_autumn[,c(2,3,4)],)
data<-data2$estimate
data.p<-data2$p.value




