library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
#setwd("/Users/ryannguyen/Desktop")

### This imports in your OTU shared file, make sure it's a csv
otu_mat <- data.frame(read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared.shared.csv",
                               header=T, 
                               sep="\t"))
otu_mat <- cbind(otu_mat[2], otu_mat[4:length(otu_mat)])

otu_mat <- t(otu_mat)

colnames(otu_mat)<- otu_mat[1,]
otu_mat <- otu_mat[-1,]
rownames(otu_mat)<-paste0("OTU", 1:nrow(otu_mat))

class(otu_mat) <- "numeric"

### This is for your taxonomy files
tax_mat <- data.frame(read.csv("stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy.csv",
                               header=T, 
                               sep="\t"))
subset_taxmat <- gsub(pattern = ";", replacement = "\t", tax_mat$Taxonomy,fixed =T)
subset_taxmat <- do.call(rbind.data.frame, strsplit(subset_taxmat, split = "\t"))


colnames(subset_taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(subset_taxmat) <- paste0("OTU", 1:nrow(subset_taxmat))
tax_mat <- subset_taxmat
tax_mat["OTU_ID"] <- NA
tax_mat["OTU_ID"] <- paste0(rownames(subset_taxmat),subset_taxmat$Genus,sep="_")

tax_mat<-as.matrix(tax_mat)
OTU <- otu_table(otu_mat,taxa_are_rows = T)
TAX <- tax_table(tax_mat)

Soil_physeq = phyloseq(OTU,TAX)

Soil_physeq <- prune_taxa(taxa_sums(Soil_physeq) > 0, Soil_physeq)
sum(taxa_sums(Soil_physeq) == 0)
#Transform to proportions
Soil_physeqP = transform_sample_counts(Soil_physeq, function(x) 100 * x/sum(x))

#Plot the top40 OTUs
top20otus = names(sort(taxa_sums(Soil_physeqP), TRUE)[1:40])
taxtab20 = cbind(tax_table(Soil_physeqP), genus40 = NA)
taxtab20[top20otus, "genus40"] <- as(tax_table(Soil_physeqP)[top20otus, "OTU_ID"], "character")
tax_table(Soil_physeqP) <- tax_table(taxtab20)

Soil_physeq_P20 = prune_taxa(top20otus, Soil_physeqP)
title = "Relative Abundance of Top40 OTUs"
pdf(file = "TopOTU's.pdf",width =20, height = 12)
plot_bar(Soil_physeq_P20, fill = "genus40", title = title) #+ labs(colour = "genus")
dev.off()
#png(file="mygraphic.png",width=400,height=350)
#plot_bar(physeq,fill="Kingdom")+ facet_wrap(~Family)
#dev.off()
#plot_bar(physeq,fill="Kingdom")
