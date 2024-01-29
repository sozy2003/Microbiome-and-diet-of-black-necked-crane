
# alpha diversity
# clean environment variables
rm(list=ls()) 

# load related packages
library(ggplot2) 
library(car)
library(performance)
library("lme4")
library(performance)

#   Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#alpha diversity
alphadiv <-  read.table("vegan.txt", header=T, row.names= 1, sep="\t")

alphadiv <-  alphadiv[rownames(design), ]

alpha = cbind(alphadiv, design)
alpha$SamplingSeason <- factor(alpha$SamplingSeason,levels = c("Spring", "Autumn"))  

p=ggplot(alpha,aes(x=SamplingSeason,y=simpson, fill=SamplingSeason))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5 )+
  geom_point(aes(x=SamplingSeason,y=simpson,),alpha=1,color= "black")+
  labs(x="SamplingSeason", y="Simpson Index") +
  scale_fill_manual()+
  theme_classic()

p=p+scale_fill_manual(values=c("#FFCC00","#33CC33" ))
p.simpson.point=p
p.simpson.point
#ggsave(paste("alpha_point_simpson",".pdf", sep=""), p.simpson.point, width =4, height =5 )

#statistical test

vegan <- alphadiv
rownames(design)==rownames(vegan)
# generate data frame for fitting
dat <- data.frame(vegan, SamplingLocationName=design$SamplingLocationName,SamplingLocation=design$SamplingLocation, SamplingSeason=design$SamplingSeason )

#lmer
mr<-lmer(simpson~SamplingSeason+(1|SamplingLocationName), data =dat)
anova(mr)#
car::Anova(mr)
check_normality(mr) # pass  
check_homogeneity(mr)# no pass    


p1<-powerTransform(mr)
m<- lmer(bcPower(simpson, p1$roundlam)~SamplingSeason+(1|SamplingLocationName), dat, REML=T)
car::Anova(m)

#shannon
# mr<-lmer(shannon~SamplingSeason+(1|SamplingLocationName), data =dat)
# anova(mr)#
# car::Anova(mr)
# 
# check_normality(mr) # no pass  
# check_homogeneity(mr)# no pass    
# 
# p1<-powerTransform(mr)
# m<- lmer(bcPower(shannon, p1$roundlam)~SamplingSeason+(1|SamplingLocationName), dat, REML=T)
# car::Anova(m)


#richness

# mr<-glmer(richness~SamplingSeason+(1|SamplingLocationName), data =dat, family = poisson(link = "log"))
# anova(mr)
# car::Anova(mr)
