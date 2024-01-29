
##The top9 family of plant diet in black-necked crane


# clean environment variables
rm(list=ls()) 

#install and load the packages
if(!require(reshape2))install.packages("reshape2")
if(!require(ggalluvial))install.packages("ggalluvial")

library("reshape2", quietly=T, warn.conflicts=F)
library("ggalluvial")
library(tidyr)
library(lme4)
library(performance)

# Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 



design$phase <- design$SamplingSeason

design$phase[design$Spring] <- "p1"  
design$phase[design$Autumn] <- "p2"

tax_count<-read.table("rbcl_otu_tax.txt",sep = "\t",header=T)


tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[6], FUN=sum) # mean

rownames(tax_count_sum) = tax_count_sum$family

tax_count_sum = tax_count_sum[,-1]
tax_count_sum = tax_count_sum[-119,]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 



mean_sort = per[(order(-rowSums(per))), ] 
colSums(mean_sort)

mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(10-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[10] = c("Low Abundance")

# write.table(mean_sort, file="Top10phylum_ProClass.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)

topN=rownames(mean_sort)


design <- design[rownames(design) %in%colnames(mean_sort),]
mean_sort <- mean_sort[,rownames(design)]


mat=mean_sort

mat_t = t(mat)

mat_t2 = merge(design[c("SamplingSeason")], mat_t, by="row.names")

mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]

geno = mat_mean$SamplingSeason
colnames(mat_mean_final) = geno

mat_mean_final = as.data.frame(mat_mean_final)

mat_mean_final$family = rownames(mat_mean_final)

data_all = as.data.frame(melt(mat_mean_final, id.vars=c("family")))

data_all$variable <- factor(data_all$variable,levels = c("Spring", "Autumn"))


data_all$value<-as.numeric(data_all$value)

# alluvium
ax<- separate(data_all,family, c("A", "family"),sep = "__")
ax[is.na(ax)]<- "Low Abundance"

ax<- ax[,-1]
data_all<- ax


data_all$family <- factor(data_all$family,levels = c("Rosaceae", "Apiaceae","Cyperaceae"
                                                     ,"Polygonaceae","Potamogetonaceae"
                                                     ,"Elaeagnaceae","Juncaceae",
                                                     "Entodontaceae","Ranunculaceae","Low Abundance"))



p = ggplot( data_all, aes(x = variable, y = value, alluvium = family, stratum=family)) +
  geom_alluvium(aes(fill = family), alpha = 0.75) +
  geom_stratum(aes(fill=family))+
  
  labs(x="SamplingSeason", y="Relative Abundance (%)")+
  theme_classic()
tax.alluvium=p

tax.alluvium
# ggsave("tax_alluvium_family_top9_rbcl.pdf", tax.alluvium, width = 5, height = 4)



mat_t2_melt <- melt(mat_t2, id.vars = "SamplingSeason")
colnames(mat_t2_melt)[2] <- "Family"
mat_t2_melt$SamplingSeason <- factor(mat_t2_melt$SamplingSeason,levels = c("Spring", "Autumn")) 
top10family_relative <- ggplot(mat_t2_melt,aes(x=SamplingSeason,y=value, fill=Family))+
  geom_boxplot()+
  geom_jitter(size=0.5, alpha=0.5)+
  facet_wrap(~Family)+theme_classic()+
  labs(x="SamplingSeason", y="Relative Abundance (%)")
top10family_relative
# ggsave("top10family_relative_boxplot_rbcl.pdf",top10family_relative, width = 10, height = 7 )

# statistics 

spp<-tax_count_sum
colSums(spp)
spp_sub <- spp[,rownames(design)]
mat_t <- t(spp_sub)
mat_t4 <-  merge(design[c("SamplingSeason", "SamplingLocationName")], mat_t, by="row.names")
rownames(mat_t4) <- mat_t4$Row.names
mat_t4 <- mat_t4[,-1]
mat_t4<-mat_t4[rownames(mat_t),]
mat_t4$total_count <- rowSums(mat_t)
table(mat_t4$SamplingLocationName)



fit <- glmer(cbind(f__Rosaceae,(total_count-f__Rosaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p1<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Apiaceae,(total_count-f__Apiaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p2<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Cyperaceae,(total_count-f__Cyperaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p3<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Polygonaceae,(total_count-f__Polygonaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p4<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Potamogetonaceae,(total_count-f__Potamogetonaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p5<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Elaeagnaceae,(total_count-f__Elaeagnaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p6<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Juncaceae,(total_count-f__Juncaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p7<-car::Anova(fit)
r2(fit)

fit <- glmer(cbind(f__Entodontaceae,(total_count-f__Entodontaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p8<-car::Anova(fit)
r2(fit)


fit <- glmer(cbind(f__Ranunculaceae,(total_count-f__Ranunculaceae))~SamplingSeason+(1|SamplingLocationName), mat_t4, family=binomial)
summary(fit)
p9<-car::Anova(fit)
r2(fit)



final<-rbind(p1,p2,p3,
             p4,p5,p6,
             p7,p8,p9)
FDR<- p.adjust(final$`Pr(>Chisq)`, method = "BH")
final<- cbind(final,FDR) 

#write.table(final,file = "final_rbcl.txt",sep = "\t",quote = F,row.names = F)
