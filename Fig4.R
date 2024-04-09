
#co-occurrence network

# clean environment variables
rm(list=ls())

library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(Biostrings)
library(phyloseq)
library(ggClusterNet)
library(tidyverse)
library(Biostrings)
library(readxl)
library(ggraph)
library(reshape2)
library(ggplot2)

##Fig4A##
Envnetplot<- paste("./16S_coi_network_notselect_one",sep = "")
dir.create(Envnetplot)


otutabcoi = read.delim("coi_otu_tax.txt", row.names=1)
otutabcoi<-otutabcoi[,-c(1:7)]
taxonomycoi = read.delim("coi_otu_tax.txt", row.names=1)
taxonomycoi<-taxonomycoi[,c(1:7)]

metadatacoi = read.delim("metadata.txt",row.names = 1)
metadatacoi<-metadatacoi[colnames(otutabcoi),]

metadatacoi$Group<-metadatacoi$SamplingSeason

#microbiota
metadata = read.delim("metadata.txt",row.names = 1)
metadata$Group<-metadata$SamplingSeason
metadata<-metadata[colnames(otutabcoi),]


otutab = read.delim("otutab.txt", row.names=1)
otutab<-otutab[,colnames(otutabcoi)]
taxonomy = read.delim("taxonomy.txt", row.names=1)

colnames(taxonomycoi)<-colnames(taxonomy)

psITS = phyloseq(sample_data(metadatacoi),
                 otu_table(as.matrix(otutabcoi), taxa_are_rows=TRUE),
                 tax_table(as.matrix(taxonomycoi)))

ps16s = phyloseq(sample_data(metadata),
                 otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
                 tax_table(as.matrix(taxonomy)))

ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16s,
                                       psITS = psITS,
                                       N16s = 500,
                                       NITS = 500,
                                       dat2.lab = "COI"
)

ps.merge

map =  phyloseq::sample_data(ps.merge)

head(map)
 map$Group = "animal"
phyloseq::sample_data(ps.merge) <- map

result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        
                        r.threshold = 0.8, 
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, 
                        # envGroup = Gru,
                        # layout = "fruchtermanreingold",
                        path = Envnetplot,
                        fill = "Phylum", 
                        size = "igraph.degree", 
                        scale = TRUE, 
                        bio = TRUE, 
                        zipi = F, 
                        step = 100, 
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)

tem <- model_maptree(cor =result[[5]],
                     method = "cluster_fast_greedy",
                     seed = 12
)
node_model = tem[[2]]
head(node_model)

p = result[[1]]
p


data = result[[2]]
plotname1 = paste(Envnetplot,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 20,height = 19)
tablename <- paste(Envnetplot,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)
tablename <- paste(Envnetplot,"/node_model_imformation",".csv",sep = "")
write.csv(node_model,tablename)

tablename <- paste(Envnetplot,"/nodeG_plot",".csv",sep = "")
write.csv(result[[4]],tablename)
tablename <- paste(Envnetplot,"/edge_plot",".csv",sep = "")
write.csv(result[[3]],tablename)
tablename <- paste(Envnetplot,"/cor_matrix",".csv",sep = "")
write.csv(result[[5]],tablename)


###Fig4B##
Envnetplot<- paste("./16S_rbcl_network_notselect_one",sep = "")
dir.create(Envnetplot)


otutabrbcl = read.delim("rbcl_otu_tax.txt", row.names=1)
otutabrbcl<-otutabrbcl[,-c(1:7)]
taxonomyrbcl = read.delim("rbcl_otu_tax.txt", row.names=1)
taxonomyrbcl<-taxonomyrbcl[,c(1:7)]

metadatarbcl = read.delim("metadata.txt",row.names = 1)
metadatarbcl<-metadatarbcl[colnames(otutabrbcl),]

metadatarbcl$Group<-metadatarbcl$SamplingSeason

metadatarbcl = read.delim("metadata.txt",row.names = 1)
metadatarbcl<-metadatarbcl[colnames(otutabrbcl),]

metadatarbcl$Group<-metadatarbcl$SamplingSeason

#microbiota
metadata = read.delim("metadata.txt",row.names = 1)
metadata$Group<-metadata$SamplingSeason
metadata<-metadata[colnames(otutabrbcl),]


otutab = read.delim("otutab.txt", row.names=1)
otutab<-otutab[,colnames(otutabrbcl)]
taxonomy = read.delim("taxonomy.txt", row.names=1)

colnames(taxonomyrbcl)<-colnames(taxonomy)

psITS = phyloseq(sample_data(metadatarbcl),
                 otu_table(as.matrix(otutabrbcl), taxa_are_rows=TRUE),
                 tax_table(as.matrix(taxonomyrbcl)))

ps16s = phyloseq(sample_data(metadata),
                 otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
                 tax_table(as.matrix(taxonomy)))



ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16s,
                                       psITS = psITS,
                                       N16s = 500,
                                       NITS = 500,
                                       dat2.lab = "RbcL"
)

ps.merge

map =  phyloseq::sample_data(ps.merge)

head(map)
map$Group = "plant"
phyloseq::sample_data(ps.merge) <- map


result <- corBionetwork(ps = ps.merge,
                        N = 0,
                        
                        r.threshold = 0.8, 
                        p.threshold = 0.05,
                        group = "Group",
                        # env = data1, 
                        # envGroup = Gru,
                        # layout = "fruchtermanreingold",
                        path = Envnetplot,
                        fill = "Phylum", 
                        size = "igraph.degree", 
                        scale = TRUE, 
                        bio = TRUE, 
                        zipi = F, 
                        step = 100, 
                        width = 12,
                        label = TRUE,
                        height = 10,
                        big = TRUE,
                        select_layout = TRUE,
                        layout_net = "model_maptree2",
                        clu_method = "cluster_fast_greedy"
                        
                        
)

tem <- model_maptree(cor =result[[5]],
                     method = "cluster_fast_greedy",
                     seed = 12
)
node_model = tem[[2]]
head(node_model)

p = result[[1]]
p

data = result[[2]]
plotname1 = paste(Envnetplot,"/network_all.pdf",sep = "")
ggsave(plotname1, p,width = 20,height = 19)
tablename <- paste(Envnetplot,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)
tablename <- paste(Envnetplot,"/node_model_imformation",".csv",sep = "")
write.csv(node_model,tablename)

tablename <- paste(Envnetplot,"/nodeG_plot",".csv",sep = "")
write.csv(result[[4]],tablename)
tablename <- paste(Envnetplot,"/edge_plot",".csv",sep = "")
write.csv(result[[3]],tablename)
tablename <- paste(Envnetplot,"/cor_matrix",".csv",sep = "")
write.csv(result[[5]],tablename)



##Fig4c##

animal<-read.csv("16S_coi_network_notselect_one/co-occurrence_Grobel_net.csv",row.names = 1)

plant<-read.csv("16S_rbcl_network_notselect_one/co-occurrence_Grobel_net.csv",row.names = 1)

data_all<-cbind(animal,plant)
data_all$net_index<-row.names(data_all)

data_all=data_all[which(rowSums(data_all[,c(1,2)]) > 0),]


data_all = as.data.frame(melt(data_all, id.vars=c("net_index")))
colnames(data_all)[2]<-"diet"


data_all2<- data_all[data_all$net_index %in% c("num.vertices","num.edges","connectance","average.degree"),]
colnames(data_all2)[2]<-"Diet"
data_all2[data_all2 == 'num.vertices'] <- 'num.nodes'
data_all2$net_index <- factor(data_all2$net_index,levels = c("num.nodes", "num.edges","connectance","average.degree"))

data_all2
p<-ggplot(data_all2,aes(x = Diet, y = value,fill=Diet))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  facet_wrap(~net_index,scales="free_y",nrow = 1)+
  scale_fill_manual(values = c("#CC99CC", "#99CC66"))+
  labs(x="Diet",y="")
p
# ggsave("allseason_four_animal-vs-plant.pdf",p, width = 7, height = 3 )









