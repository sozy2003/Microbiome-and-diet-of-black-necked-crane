
#β diversity of plant diet

# clean environment variables
rm(list=ls()) 

# load related packages
library("ggplot2") 
library(ggsignif)
library("vegan")
library(tidyr)

library(performance)
library("lme4")
library(car)
library(sqldf)

# Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t")


otu_Flattening<-read.table("otu_rare_rbcl.txt", header=T, row.names= 1, sep="\t") 
otu_Flattening<-t(otu_Flattening)

betadiv_rbcl <- vegdist(otu_Flattening,method = "bray")
betadiv_rbcl<- as.matrix(betadiv_rbcl)
design_rbcl <-  design[row.names(design)%in%row.names(betadiv_rbcl),] 
betadiv_rbcl =betadiv_rbcl[rownames(design_rbcl), rownames(design_rbcl)] 

betadiv<-betadiv_rbcl

betadiv = betadiv[rownames(design_rbcl), rownames(design_rbcl)]

betadis <- as.matrix(betadiv)

betadis[upper.tri(betadis,diag = T)] <- NA


SamplingSeason <- unique(sort(design_rbcl$SamplingSeason))

p=list() 
for (i in SamplingSeason){
  p[[i]] <- betadis[design_rbcl$SamplingSeason==i, design_rbcl$SamplingSeason==i]
  nrow <- dim(p[[i]])[1]
  ncol <- dim(p[[i]])[2]
  row <- rep(rownames(p[[i]]), ncol)
  col <- rep(colnames(p[[i]]),each=nrow)
  p[[i]] <- data.frame(row, col, value=as.numeric(p[[i]]))
  p[[i]]$SamplingSeason <- rep(i,nrow(p[[i]]))
  
}
p



frame <- rbind(p[[1]],p[[2]])


frame <- na.omit(frame)
frame$SamplingSeason <- factor(frame$SamplingSeason,levels = c("Spring", "Autumn")) 

p=ggplot(frame,aes(x=as.factor(SamplingSeason),y=value,fill=SamplingSeason))+
  geom_boxplot()+
  labs(x="SamplingSeason", y="bray_curtis") + theme_classic()
p=p+geom_jitter(position=position_jitter(0.17), size=0.001, alpha=0.3)
  # geom_smooth(method = "lm", formula = y ~ poly(x,3))
p=p+scale_fill_manual(values=c("#FFCC00","#33CC33" ))+
  geom_signif(comparisons=list(c("Spring", "Autumn")), map_signif_level=TRUE,color="black")
p.boxplot=p
p.boxplot
# ggsave(paste("beta_box_bray_curtis_rbcl",".pdf", sep=""), p.boxplot, width = 2.5, height =2.5 )





# statistics 
df<-design_rbcl[,c("SamplingLocation","SamplingLocationName")]
df$row<-row.names(df)
df2<-merge(df,frame,by="row")


colnames(df)<-c("SamplingLocation2","SamplingLocationName2","col")
df2<-merge(df,df2,by="col")

zda<-df2[which(df2$SamplingLocationName2 == df2$SamplingLocationName),]

df2<-zda






m2 <- lmer(value~SamplingSeason+(1|SamplingLocationName), df2, REML=F)
car::Anova(m2)

check_normality(m2) #no pass
check_homogeneity(m2)#no pass


p1<-powerTransform(m2)
m<- lmer(bcPower(value, p1$roundlam)~SamplingSeason+(1|SamplingLocationName), df2, REML=T)
car::Anova(m)















