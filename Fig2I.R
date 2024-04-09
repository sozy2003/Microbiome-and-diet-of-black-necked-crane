# plant diet CAP beased on bray_curtis

# clean environment variables
rm(list=ls()) 
library(vegan)
library(phyloseq)
library(ggplot2)
library(plyr)

design <- read.table("metadata.txt",row.names= 1,  header=T, sep="\t")  
all_bDAT <- read.delim("otu_rare_rbcl.txt", row.names= 1,  header=T, sep="\t")
design<-design[colnames(all_bDAT),]
bDAT <- all_bDAT[, rownames(design) ]

bTAX <- read.table( "taxonomy_rbcl.txt" , row.names=1, sep="\t", header=T)


bTAX <- bTAX[rownames(bDAT),]

bTAX$OTU_ID <- rownames(bTAX)

bDAT_field <- bDAT[, rownames(design) ]


bDAT_field_rare<-bDAT_field
bDAT_field_rare <- bDAT_field_rare[rowSums(bDAT_field_rare) >= 1,]   


field_phy <- phyloseq(sample_data(design), 
                      otu_table(bDAT_field_rare, taxa_are_rows=T), 
                      tax_table(as.matrix(bTAX[rownames(bDAT_field_rare),])) )


head(design)
field_CAP_bray <- ordinate(field_phy, method="CAP", ~SamplingSeason, distance="bray")


variability_table <- function(cca){
  
  chi <- c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi / chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained") 
  return(variability_table)
  
}

field_CAP_bray_var_tbl <- variability_table(field_CAP_bray)

# ANOVA
set.seed(000)
field_CAP_bray_anova <- anova.cca(field_CAP_bray, permutations=999)

# plot
title <- paste("[~SamplingSeason: ",
               format(field_CAP_bray_var_tbl["constrained", "proportion"] * 100, digits=2, nsmall=1),
               "% of variance, P = ",
               format(field_CAP_bray_anova[1, 4], digits=2),
               "]", sep = "")
 pdf("Field_CAP_season_rbcl.pdf" ,width=5, height=3,)
p2 <- plot_ordination(field_phy, field_CAP_bray, color ="SamplingSeason" )
p2 <- p2 + geom_point(size=3)
p2 <- p2 + scale_color_manual(values=c("#33CC33","#FFCC00" ))
 p2 <- p2 + ggtitle(title)
p2 <- p2 + theme(plot.title = element_text(size=9, face = "bold"))
p2

dev.off()
print(p2)

# PERMANOVA
bDAT_field_rare_dist <- vegdist(t(bDAT_field_rare), method="bray")
set.seed(000)
bDAT_field_rare_dist_paov <- adonis2(bDAT_field_rare_dist ~SamplingSeason, data=design, permutations=999)
bDAT_field_rare_dist_paov

