
# α diversity of animal diet

# clean environment variables
rm(list=ls()) 

# load related packages
if(!require(ggplot2))install.packages("ggplot2")
library("ggplot2") 
library(ggsignif)
library(tidyr)
library(vegan)
library(car)
library(performance)
library("lme4")
#   Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#alpha diversity

otu_Flattening<- read.table("otu_rare_coi.txt", header=T, row.names= 1, sep="\t")
otu_Flattening<-t(otu_Flattening)
#α diversity
shannon <- diversity(otu_Flattening, index="shannon", MARGIN = 1)
simpson <- diversity(otu_Flattening, index="simpson", MARGIN = 1)
richness <- rowSums( otu_Flattening != 0)
alphadiv_coi_fla <- data.frame(richness=richness, shannon=shannon, simpson=simpson)
alphadiv <-alphadiv_coi_fla

design<-design[rownames(design)%in% rownames (alphadiv),]
alphadiv <-  alphadiv[rownames(design), ]
# write.table(alphadiv,"alphadiv_coi.txt",row.names= T, sep="\t")

alpha = cbind(alphadiv, design)
alpha$SamplingSeason <- factor(alpha$SamplingSeason,levels = c("Spring", "Autumn")) 

# simpson
p=ggplot(alpha,aes(x=SamplingSeason,y=simpson,fill=SamplingSeason))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5)+
  geom_point(aes(x=SamplingSeason,y=simpson,),alpha=1)+
  labs(x="SamplingSeason", y="Simpson Index") +
  scale_fill_manual()+
  theme_classic()

p=p+scale_fill_manual(values=c("#FFCC00","#33CC33" ))
   # geom_signif(comparisons=list(c("Spring", "Autumn")), map_signif_level=TRUE,color="black")
p 

# ggsave(paste("alpha_box_simpson_coi",".pdf", sep=""), p, width = 3, height =3 )

# richness
# p=ggplot(alpha,aes(x=SamplingSeason,y=richness,color=SamplingSeason))+
#   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent")+
#   geom_point(aes(x=SamplingSeason,y=richness,),alpha=1)+
#   labs(x="SamplingSeason", y="Richness Index") +
#   scale_fill_manual()+
#   theme_classic()
# 
# p=p+scale_color_manual(values=c("#FFCC00","#33CC33" ))
# p

#ggsave(paste("alpha_point_richness",".pdf", sep=""), p, width =3.5, height =4 )

# p=ggplot(alpha,aes(x=SamplingSeason,y=shannon,color=SamplingSeason))+
#   geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent")+
#   geom_point(aes(x=SamplingSeason,y=shannon,),alpha=1)+
#   labs(x="SamplingSeason", y="Shannon Index") +
#   scale_fill_manual()+
#   theme_classic()
# 
# p=p+scale_color_manual(values=c("#FFCC00","#33CC33" ))
# p.shannon.point=p
# p.shannon.point
#ggsave(paste("alpha_point_shannon",".pdf", sep=""), p.shannon.point, width =4, height =5 )


# statistics 
#richness
# mr<-glmer(richness~SamplingSeason+(1|SamplingLocationName), data =alpha, family = poisson(link = "log"))
# anova(mr)
# car::Anova(mr)


#shannon
# mr<-lmer(shannon+1~SamplingSeason+(1|SamplingLocationName), data =alpha)
# anova(mr)
# car::Anova(mr)
# 
# check_normality(mr) # no pass  
# check_homogeneity(mr)# no pass    
# p1<-powerTransform(mr)
# m<- lmer(bcPower(shannon+1, p1$roundlam)~SamplingSeason+(1|SamplingLocationName), alpha, REML=T)
# car::Anova(m)


#simpson
mr<-lmer(simpson+1~SamplingSeason+(1|SamplingLocationName), data =alpha)
anova(mr)#
car::Anova(mr)


check_normality(mr) # no pass 
check_homogeneity(mr)#  pass  


p1<-powerTransform(mr)
m<- lmer(bcPower(simpson+1, p1$roundlam)~SamplingSeason+(1|SamplingLocationName), alpha, REML=T)
car::Anova(m)

