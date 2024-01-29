
#plant diet NMDS beased on bray_curtis
# clean enviroment object
rm(list=ls()) 

library(vegan)
library(ggplot2)
# Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t")

design$phase <- design$SamplingSeason

design$phase[design$Spring] <- "p1"
design$phase[design$Autumn] <- "p2"
# bray_curtis

library(tidyr)
otu_Flattening<-read.table("otu_rare_rbcl.txt", header=T, row.names= 1, sep="\t") 
otu_Flattening<-t(otu_Flattening)

betadiv_rbcl <- vegdist(otu_Flattening,method = "bray")
betadiv_rbcl<- as.matrix(betadiv_rbcl)
design_rbcl <-  design[row.names(design)%in%row.names(betadiv_rbcl),] 




betadiv_rbcl = betadiv_rbcl[rownames(design_rbcl), rownames(design_rbcl)]


set.seed(000)
spe.nmds <- metaMDS(betadiv_rbcl,distance="bray")
spe.nmds
df_nmds_stress <-spe.nmds$stress
df_nmds_stress


df_points <- as.data.frame(spe.nmds$points)

df_points$samples <- row.names(df_points)

names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)

library("dplyr")
df_points<-df_points %>%
  mutate(season = case_when(
    grepl("^A\\d", samples) ~ "Spring",
    grepl("^B\\d",samples) ~ "Autumn",))

df <-df_points

df$season <- factor(df$season,levels = c("Spring", "Autumn"))  
color=c("#FFCC00","#33CC33")
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+
  theme_bw()+
  geom_point(aes(color = season), shape = 19, size=3)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_text(aes(label=samples, y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0, color = season),size=3.5, show.legend = F)+
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=season),
               alpha=0.2)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())+
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))
p1
#ggsave("NMDS_BC_rbcl_2023.7.28.pdf", p1, width = 5, height = 4)
set.seed(000)
anosim=anosim(betadiv_rbcl,design_rbcl$SamplingSeason, permutations=999)
anosim
summary(anosim)
#figS6
anosim$class.vec<-factor(anosim$class.vec,levels = c("Between","Spring","Autumn"))

anosim2<-plot(anosim,col=c("grey","yellow","green"))
 
 pdf("NMDS_anosim_BC_rbcl.pdf")

 plot(anosim,col=c("grey","yellow","green"))

 dev.off()
 