
# The relationship between microbial and diet alpha diversity. 

# clean environment variables
rm(list=ls()) 

#16s

# load related packages
if(!require(ggpmisc))install.packages("ggpmisc")
library("ggplot2") 
library("ggpmisc")
library(vegan) 
library(tidyr)



#   Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#alpha diversity
alphadiv <-  read.table("vegan.txt", header=T, row.names= 1, sep="\t")
alphadiv <-  alphadiv[rownames(design), ]

alpha = cbind(alphadiv, design)
alpha <- alpha[,c(1:6,10)]

#rbcl
design_rbcl <-  design[-c(1,3),] 

#alpha diversity





otu_Flattening<-read.table("otu_rare_rbcl.txt", header=T, row.names= 1, sep="\t") 
otu_Flattening<-t(otu_Flattening)

shannon <- diversity(otu_Flattening, index="shannon", MARGIN = 1)
shannon
simpson <- diversity(otu_Flattening, index="simpson", MARGIN = 1)
simpson

richness <- rowSums( otu_Flattening != 0)
richness


alphadiv_rbcl <- data.frame(richness=richness, shannon=shannon, simpson=simpson)

alphadiv_rbcl <-  alphadiv_rbcl[rownames(design_rbcl), ]


alpha_rbcl = cbind(alphadiv_rbcl, design_rbcl)
alpha_rbcl <- alpha_rbcl[,c(1:3,7)]
head(alpha_rbcl)
colnames(alpha_rbcl) <- c("richness_rbcl","shannon_rbcl","simpson_rbcl","SamplingSeason_rbcl")

d <- alpha
names <- rownames(d)
rownames(d) <- NULL
alpha <- cbind(names,d)

d <- alpha_rbcl
names <- rownames(d)
rownames(d) <- NULL
alpha_rbcl <- cbind(names,d)

frame_rbcl<-merge(alpha_rbcl,alpha,by="names")
frame_rbcl$SamplingSeason <- factor(frame_rbcl$SamplingSeason,levels = c("Spring", "Autumn"))




# plot
p4 <- ggplot(frame_rbcl,aes(x = simpson_rbcl, y = simpson))+
  
  geom_point(aes(color=SamplingSeason))+
  
  geom_smooth(method = "lm")+
  scale_y_continuous(name = "Simpson microbiome")+
  scale_x_continuous(name = "Simpson(rbcl)")+
  geom_smooth(frame_rbcl[frame_rbcl$SamplingSeason=="Spring",], mapping=aes(x = simpson_rbcl, y = simpson,color=SamplingSeason),method = "lm")+
geom_smooth(frame_rbcl[frame_rbcl$SamplingSeason=="Autumn",], mapping=aes(x = simpson_rbcl, y = simpson,color=SamplingSeason),method = "lm")

p4<-p4+scale_color_manual(values=c("#FFCC00","#33CC33" ))+
  annotate("text",label="pearson:r==0.2823~p<0.05",family="Times",parse=T,x=0.3,y=0.25,color="blue",size=4)
  
p4

# ggsave("pearson_rbcl_richness_diet_micro.pdf", p4, width =6, height = 4)

#coi


#alpha diversity


otu_Flattening<-read.table("otu_rare_coi.txt", header=T, row.names= 1, sep="\t") 
otu_Flattening<-t(otu_Flattening)

shannon <- diversity(otu_Flattening, index="shannon", MARGIN = 1)
shannon
simpson <- diversity(otu_Flattening, index="simpson", MARGIN = 1)
simpson

richness <- rowSums( otu_Flattening != 0)
richness


alphadiv_coi <- data.frame(richness=richness, shannon=shannon, simpson=simpson)


design_coi <-  design[row.names(design)%in%row.names(alphadiv_coi),] 


alphadiv_coi <-  alphadiv_coi[rownames(design_coi), ]


alpha_coi = cbind(alphadiv_coi, design_coi)
alpha_coi <- alpha_coi[,c(1:3,7)]
head(alpha_coi)
colnames(alpha_coi) <- c("richness_coi","shannon_coi",
                         "simpson_coi","SamplingSeason_coi")


d <- alpha_coi
names <- rownames(d)
rownames(d) <- NULL
alpha_coi <- cbind(names,d)

frame_coi<-merge(alpha_coi,alpha,by="names")
frame_coi$SamplingSeason <- factor(frame_coi$SamplingSeason,levels = c("Spring", "Autumn"))


p8 <- ggplot(frame_coi,aes(x = simpson_coi, y =simpson))+
  
  geom_point(aes(color=SamplingSeason))+
  geom_smooth(frame_coi[frame_coi$SamplingSeason=="Spring",], mapping=aes(x = simpson_coi, y = simpson,color=SamplingSeason),method = "lm")+
  geom_smooth(frame_coi[frame_coi$SamplingSeason=="Autumn",], mapping=aes(x = simpson_coi, y = simpson,color=SamplingSeason),method = "lm")+
  geom_smooth(method = "lm")+
  scale_y_continuous(name = "Simpson microbiome")+
  scale_x_continuous(name = "Simpson(coi)")


p8<-p8+scale_color_manual(values=c("#FFCC00","#33CC33" ))+
  # stat_poly_eq(aes(label=paste(..eq.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  annotate("text",label="pearson:r==0.0066~p>0.05",family="Times",parse=T,x=0.3,y=0.25,color="blue",size=4)
  

p8
# ggsave("pearson_coi_richness_diet_micro.pdf", p8, width =6, height = 4)



library(patchwork)

p8<- p8+ ggtitle('C')
p4<- p4+ ggtitle('D')
p9<-p8+p4
p9
#ggsave("pearson_diet_micro.pdf", p9, width =10, height = 4)



# pearson coefficient 
cor.test(frame_rbcl$simpson_rbcl, frame_rbcl$simpson, method = "pearson")
cor.test(frame_coi$simpson_coi, frame_coi$simpson, method = "pearson")

#rbcl
#spring
cor.test(frame_rbcl[c(1:29),]$simpson_rbcl, frame_rbcl[c(1:29),]$simpson, method = "pearson")

#autumn
cor.test(frame_rbcl[c(30:58),]$simpson_rbcl, frame_rbcl[c(30:58),]$simpson, method = "pearson")

#coi
#spring
cor.test(frame_coi[c(1:28),]$simpson_coi, frame_rbcl[c(1:28),]$simpson, method = "pearson")

#autumn
cor.test(frame_coi[c(29:56),]$simpson_coi, frame_rbcl[c(29:56),]$simpson, method = "pearson")



#animal and plant
cor.test(frame_rbcl$richness_rbcl, frame_rbcl$richness, method = "pearson")
cor.test(frame_coi$richness_coi, frame_coi$richness, method = "pearson")


#spearman

#animal##
cor.test(frame_coi$simpson_coi, frame_coi$simpson, method = "spearman")

#spring
cor.test(frame_coi[c(1:28),]$simpson_coi, frame_coi[c(1:28),]$simpson, method = "spearman")

#autumn
cor.test(frame_coi[c(29:56),]$simpson_coi, frame_coi[c(29:56),]$simpson, method = "spearman")


#plant###
cor.test(frame_rbcl$simpson_rbcl, frame_rbcl$simpson, method = "spearman")

#spring
cor.test(frame_rbcl[c(1:29),]$simpson_rbcl, frame_rbcl[c(1:29),]$simpson, method = "spearman")

#autumn
cor.test(frame_rbcl[c(30:58),]$simpson_rbcl, frame_rbcl[c(30:58),]$simpson, method = "spearman")





#fig s9



#delect the outlier point 

frame_coi_sec<-frame_coi[-54,]

ps1 <- ggplot(frame_coi_sec,aes(x = richness_coi, y =richness))+
  
  geom_point(aes(color=SamplingSeason))+
  geom_smooth(frame_coi[frame_coi$SamplingSeason=="Spring",], mapping=aes(x = richness_coi, y = richness,color=SamplingSeason),method = "lm")+
  geom_smooth(frame_coi[frame_coi$SamplingSeason=="Autumn",], mapping=aes(x = richness_coi, y = richness,color=SamplingSeason),method = "lm")+
  geom_smooth(method = "lm")+
  scale_y_continuous(name = "Richness microbiome")+
  scale_x_continuous(name = "Richness(coi)")


ps1<-ps1+scale_color_manual(values=c("#FFCC00","#33CC33" ))+
  # stat_poly_eq(aes(label=paste(..eq.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  annotate("text",label="pearson:r==-0.3145~p<0.05",family="Times",parse=T,x=17,y=0.25,color="blue",size=4)


ps1

frame_rbcl_sec<-frame_rbcl[-7,]

ps2 <- ggplot(frame_rbcl_sec,aes(x = richness_rbcl, y =richness))+
  
  geom_point(aes(color=SamplingSeason))+
  geom_smooth(frame_rbcl[frame_rbcl$SamplingSeason=="Spring",], mapping=aes(x =richness_rbcl, y =richness,color=SamplingSeason),method = "lm")+
  geom_smooth(frame_rbcl[frame_rbcl$SamplingSeason=="Autumn",], mapping=aes(x =richness_rbcl, y =richness,color=SamplingSeason),method = "lm")+
  geom_smooth(method = "lm")+
  scale_y_continuous(name = "Richness microbiome")+
  scale_x_continuous(name = "Richness(rbcl)")


ps2<-ps2+scale_color_manual(values=c("#FFCC00","#33CC33" ))+
  # stat_poly_eq(aes(label=paste(..eq.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)+
  annotate("text",label="pearson:r==-2931~p<0.05",family="Times",parse=T,x=25,y=0.25,color="blue",size=4)


ps2

figs<-ps1+ps2
figs

#ggsave("pearson_diet_micro_richness.pdf", figs, width =10, height = 4)


#spearman

#animal##
cor.test(frame_coi_sec$richness_coi, frame_coi_sec$richness, method = "spearman")

#spring
cor.test(frame_coi_sec[c(1:28),]$richness_coi, frame_coi_sec[c(1:28),]$richness, method = "spearman")

#autumn
cor.test(frame_coi_sec[c(29:55),]$richness_coi, frame_coi_sec[c(29:55),]$richness, method = "spearman")


#plant###
cor.test(frame_rbcl_sec$richness_rbcl, frame_rbcl_sec$richness, method = "spearman")

#spring
cor.test(frame_rbcl_sec[c(1:28),]$richness_rbcl, frame_rbcl_sec[c(1:28),]$richness, method = "spearman")

#autumn
cor.test(frame_rbcl_sec[c(30:57),]$richness_rbcl, frame_rbcl_sec[c(30:57),]$richness, method = "spearman")








