
# α diversity of plant diet

# clean environment variables
rm(list=ls()) 

# load related packages
library("ggplot2") 
library(tidyr)
library(vegan)
library(car)
library(performance)
library("lme4")
#   Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#alpha diversity

otu_Flattening<- read.table("otu_rare_rbcl.txt", header=T, row.names= 1, sep="\t")
otu_Flattening<-t(otu_Flattening)

shannon <- diversity(otu_Flattening, index="shannon", MARGIN = 1)
simpson <- diversity(otu_Flattening, index="simpson", MARGIN = 1)
richness <- rowSums( otu_Flattening != 0)
alphadiv_coi_fla <- data.frame(richness=richness, shannon=shannon, simpson=simpson)
alphadiv <-alphadiv_coi_fla

design<-design[rownames(design)%in% rownames (alphadiv),]
alphadiv <-  alphadiv[rownames(design), ]

# write.table(alphadiv,"alphadiv_rbcl.txt",row.names= T, sep="\t")

alpha = cbind(alphadiv, design)
alpha$SamplingSeason <- factor(alpha$SamplingSeason,levels = c("Spring", "Autumn"))  


# plot 
p=ggplot(alpha,aes(x=SamplingSeason,y=simpson,fill=SamplingSeason))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5)+
  geom_point(aes(x=SamplingSeason,y=simpson,),alpha=1)+
  labs(x="SamplingSeason", y="Simpson Index") +
  scale_fill_manual()+
  theme_classic()

p=p+scale_fill_manual(values=c("#FFCC00","#33CC33" ))
p

# ggsave(paste("alpha_box_simpson_rbcl",".pdf", sep=""), p, width = 3, height =3 )


# statistics 
#richness
# mr<-glmer(richness~SamplingSeason+(1|SamplingLocationName), data =alpha, family = poisson(link = "log"))
# anova(mr)
# car::Anova(mr)
# mr<-lmer(shannon+1~SamplingSeason+(1|SamplingLocationName), data =alpha)
# anova(mr)
# car::Anova(mr)
# 
# check_normality(mr) # no pass  
# check_homogeneity(mr)#  pass    
# p1<-powerTransform(mr)
# m<- lmer(bcPower(shannon+1, p1$roundlam)~SamplingSeason+(1|SamplingLocationName), alpha, REML=T)
# car::Anova(m)

#simpson
mr<-lmer(simpson+1~SamplingSeason+(1|SamplingLocationName), data =alpha)
anova(mr)#
car::Anova(mr)


check_normality(mr) # pass 
check_homogeneity(mr)#  pass  

