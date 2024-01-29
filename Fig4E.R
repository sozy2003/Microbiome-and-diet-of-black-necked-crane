
# clean environment variables
rm(list=ls()) 

#Procrustes analysis
library("vegan")
library("ggplot2")
library("cowplot")

#bray-curtis similarity of diet and microbiome, Mantel test 
metadata = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

df1 =  read.table("otu_rare_coi.txt", header=T, row.names= 1, sep="\t") 

df2 = read.table("otutab_rare.txt", header=T, row.names= 1, sep="\t")
df1=t(df1)
df2=t(df2)

dist.abund <- vegdist(df1, method = "bray")
dist.abund <- as.dist(dist.abund)
mdist.abund = vegdist(df2, method = "bray")
mdist.abund <- as.dist(mdist.abund)

dpcoa <- as.data.frame(cmdscale(dist.abund)) 

mpcoa <- as.data.frame(cmdscale(mdist.abund))

mpcoa <-  mpcoa[row.names(mpcoa)%in%row.names(dpcoa),]

#procrustes analysis
pro <- procrustes(X = dpcoa, Y = mpcoa, scale = TRUE,symmetric = TRUE)


set.seed(000)
pro_test <- protest(dpcoa,mpcoa,perm=9999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Diet"

seasons=metadata[,c(4,12)]
seasons=seasons[row.names(seasons)%in%row.names(dpcoa),]
beta_pro=cbind(beta_pro,seasons)

trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome"
seasons=metadata[,c(4,12)]
seasons=seasons[row.names(seasons)%in%row.names(dpcoa),]

trans_pro=cbind(trans_pro,seasons)
colnames(trans_pro) <- colnames(beta_pro)
pval <- signif(pro_test$signif, 1)
plot <- rbind(beta_pro, trans_pro)

plot$SamplingSeason <- factor(plot$SamplingSeason,levels = c("Spring", "Autumn"))
grass_food_micro <- ggplot(plot) +
  geom_point(size = 4, alpha=0.75, aes(x = V1, y = V2, color =SamplingSeason,shape=type))+ 
  scale_color_manual(values=c("#FFCC00","#33CC33" )) +
  theme_classic() +
  scale_x_continuous()+
  scale_y_continuous()+
  geom_line(aes(x= V1, y=V2, group=UserName), col = "darkgrey", alpha = 0.6,linewidth=0.2) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10,colour="black"),
        legend.position = 'bottom',
        axis.text = element_text(size=10,colour="black"),
        axis.title = element_text(size=13,colour="black"),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) +
  annotate("text", x = 0.07, y = -0.13, label = paste0("P=",pval), size = 4) 

grass_food_micro_leg <- get_legend(grass_food_micro) 

p<-grass_food_micro + theme(legend.position = "right")
p
#ggsave("coi_microbiome_Procrustes2023.8.30.pdf", p, width = 5, height = 4)
