


#The relationship between top 10 diet and whole microbiota matrix. 

rm(list = ls()) 
#devtools::install_local("D:/R/R-4.2.1/library/ggcor", force = TRUE)

library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
library(tidyr)


#coi

##########################################################################################

tax_count<-read.table("coi_otu_tax.txt",sep = "\t",header=T)


tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[6], FUN=sum) # mean

ax<- separate(tax_count_sum,family, c("A", "family"),sep = "__")
ax<- ax[,-1]
tax_count_sum<- ax


rownames(tax_count_sum) = tax_count_sum$family

tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 

mean_sort = per[(order(-rowSums(per))), ]
colSums(mean_sort)

mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(10-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[10] = c("Low Abundance")

rownames(mean_sort)[9] = c("BOLD:ACY5875")

mean_sort <-data.frame(t(mean_sort))

#microbiome OTU table 

df <- read.table("otutab_rare.txt",sep="\t",header = T,row.names = 1,check.names = F)
df <-df[,colnames(df)%in%colnames(tax_count_sum)]
df <-data.frame(t(df))



df_mantel <- mantel_test(df, mean_sort, mantel.fun = 'mantel',
                         spec.dist.method = 'bray', 
                         env.dist.method = 'euclidean',spec.select = list(spring = 1:28,
                                                                          autumn = 29:56))
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
                    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


df_mantel$spec <- factor(df_mantel$spec,levels = c("spring", "autumn")) 
p1<-quickcor(mean_sort ,method = "spearman",exact=FALSE, type = "upper", cor.test = T, cluster.type = "all") +
  geom_square() +
  # geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+
  scale_fill_gradient2( high = '#FF3366', mid = 'white',low = '#6666CC') + 
  anno_link(df_mantel, aes(linetype = df_p,
                           size = df_r,color=spec))+
  scale_linetype()+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+
  scale_color_manual(values = c("#CC9900","#339900"))+
  guides(fill = guide_colorbar(title = "correlation", order = 1),
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "season", order = 3),
         linetype = guide_legend(title = "Mantel's p", order = 3))
p1

#ggsave("AAAmantel_coi_top10_NCBI_BOLD.pdf", p1, width = 8, height = 6)



#rbcl#######################################################
tax_count<-read.table("rbcl_otu_tax.txt",sep = "\t",header=T)
# group by family
tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[6], FUN=sum) # mean

ax<- separate(tax_count_sum,family, c("A", "family"),sep = "__")
ax<- ax[,-1]
tax_count_sum<- ax


rownames(tax_count_sum) = tax_count_sum$family

tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 

mean_sort = per[(order(-rowSums(per))), ] 
colSums(mean_sort)

mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(10-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[10] = c("Low Abundance")

mean_sort <-data.frame(t(mean_sort))

df <- read.table("otutab_rare.txt",sep="\t",header = T,row.names = 1,check.names = F)
df <-data.frame(t(df))
df <-df[-c(6,26),]


df_mantel <- mantel_test(df, mean_sort, mantel.fun = 'mantel',
                         spec.dist.method = 'bray', 
                         env.dist.method = 'euclidean',spec.select = list(spring = 1:28,
                                                                          autumn = 29:58))
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
                    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


df_mantel$spec <- factor(df_mantel$spec,levels = c("spring", "autumn")) 
p2<-quickcor(mean_sort ,method = "spearman",exact=FALSE, type = "upper", cor.test = T, cluster.type = "all") +
  geom_square() +
  # geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+
  scale_fill_gradient2(high = '#FF3366', mid = 'white',low = '#6666CC') + 
  anno_link(df_mantel, aes(linetype = df_p,
                           size = df_r,color=spec))+
  scale_linetype()+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+
  scale_color_manual(values = c("#CC9900","#339900"))+
  guides(fill = guide_colorbar(title = "correlation", order = 1),
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "season", order = 3),
         linetype = guide_legend(title = "Mantel's p", order = 3))
p2
#ggsave("mantel_rbcl_top10.pdf", p2, width = 8, height = 6)

#install.packages("devtools")
library(devtools)
#install_github("thomasp85/patchwork")
library(ggplot2)
library(patchwork)
#p<- p+ ggtitle('anml')
p1<- p1+ ggtitle('coi')
p2<- p2+ ggtitle('rbcl')
p3<-p1+p2
p3
#ggsave("mantel_anml_coi_rbcl_top10.pdf", p3, width = 10, height = 5)

