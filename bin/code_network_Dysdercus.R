#### Scripts used in the article: ####
# Transgenic-wild-cotton-and-its-unintended-effects-on-the-gut-microbiota-of-genus-Dysdercus

# Authors: Javier Pérez-López ▲‡, Valeria Alavez ▲, Juan Fornoni», René Cerritos•, Ana Wegier Φ ‡.

# ▲ Posgrado en Ciencias Biológicas, Instituto de Biología, Universidad Nacional Autónoma de México, Ciudad de México, México.
# » Instituto de Ecología, Universidad Nacional Autónoma de México, AP 70-275, Ciudad Universitaria, Coyoacán, 04510 CDMX, México.
# • División de Investigación, Facultad de Medicina, Universidad Nacional Autónoma de México, Av. Universidad 3000, circuito escolar s/n, 04510, Ciudad de México, México.
# Φ Laboratorio de Genética de la Conservación, Jardín Botánico, Instituto de Biología, Universidad Nacional Autónoma de México, Ciudad de México, México.

rm(list = ls())

# libraries 
library(microbiome)
library(ggpubr)
library(qiime2R)
library(NetCoMi)
library(igraph)
library(tidyverse)
library(microbiomeutilities)
library(dplyr)
library(patchwork)

setwd("~/Escritorio/NCBI_formato/bin/")

physeq<-qza_to_phyloseq(
  features="../files/2_table.qza",
  tree="../files/7_rooted-tree.qza",
  "../files/6_taxonomy.qza",
  metadata = "../files/Metadata_manual.csv")

physeq_filtro<-subset_samples(physeq, Especie %in% c("concinnus", "obliquus"))

physeq_hembras_neg <-subset_samples(physeq_filtro, sex == "hembra") %>%
  subset_samples(PCR_result=="Negativo")

physeq_hembras_cry <-subset_samples(physeq_filtro, sex == "hembra") %>%
  subset_samples(PCR_result=="cry1")

physeq_macho_neg <-subset_samples(physeq_filtro, sex == "macho") %>%
  subset_samples(PCR_result=="Negativo")

physeq_macho_cry <-subset_samples(physeq_filtro, sex == "macho") %>%
  subset_samples(PCR_result=="cry1")

physeq_hembras <-subset_samples(physeq_filtro, sex == "hembra") %>% 
  subset_samples(Especie %in% c("concinnus", "obliquus"))

physeq_machos <-subset_samples(physeq_filtro, sex != "hembra") %>% 
  subset_samples(Especie %in% c("concinnus", "obliquus"))

#### FIGURE 1 #####
########## Composition of gut microbiota of males and females of Dysdercus collected on wild cotton with and without transgenes #######################################
erie_phylum <- physeq_filtro %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                               # Replacing values

# Set colors for plotting
phylum_colors <- c("#d4b8ff",
                   "#b3de69",
                   "#009bb9",
                   "#f8785c",
                   "#02ceba",
                   "orange")

#### Plot Phyla level ####
levels(erie_phylum$sex) <- c("Female", "Male")
levels(erie_phylum$PCR_result) <- c("cry1ab/ac", "no-cry1ab/ac")
levels(erie_phylum$Phylum)

p1<- ggplot(erie_phylum,
            aes(x =PCR_result, y = Abundance,
                fill = Phylum)) + 
  facet_grid(sex~., scales = "free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") + xlab("")+
  ggtitle("A)") +
  theme_pubclean()+ theme(legend.position = "right")+ 
  coord_flip()+ theme(
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 8))+                               # Change margins of ggplot2 plot
  theme(plot.margin = unit(c(2.5, 2.5, 0, 2.5), "cm"))
p1

### Plot Class level #####
erie_phylum <- physeq_filtro %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)   

# Plot 
names(erie_phylum)
unique(erie_phylum$Class)
levels(erie_phylum$sex) <- c("Female", "Male")
levels(erie_phylum$PCR_result) <- c("cry1ab/ac", "no-cry1ab/ac")

p2<- ggplot(erie_phylum %>% filter(Class != "Chloroplast"), 
            aes(x =PCR_result, y = Abundance, fill = Class)) + 
  facet_grid(sex~., scales = "free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5')) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Class > 2%) \n") + xlab("")+
  ggtitle("B)")+theme_pubclean()+ 
  theme(legend.position = "right") + 
  coord_flip()+ theme(
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 8))+                               # Change margins of ggplot2 plot
  theme(plot.margin = unit(c(0, 2.5, 2.5, 2.5), "cm"))


p1 / p2

### Analysis of gut microbiota of females of Dysdercus collected on wild cotton with and without transgenes ####
### FIGURE 2A: Network female without cry1ab/ac ####
amgut_genus_neg <- phyloseq::tax_glom(physeq_hembras_neg, taxrank = "Genus")
taxtab <- amgut_genus_neg@tax_table@.Data
amgut_genus_neg@tax_table@.Data <- taxtab
rownames(amgut_genus_neg@otu_table@.Data) <- taxtab[, "Genus"]

# Network construct female whitout cry1ab/ac 
net_hembra_neg <- netConstruct(amgut_genus_neg, 
                               filtTax = "numbSamp",
                               filtTaxPar = list(numbSamp = 1),
                               filtSamp = "totalReads",
                               filtSampPar = list(numbTaxa = 50),
                               measure = "pearson",
                               zeroMethod = "multRepl",
                               normMethod = "clr", 
                               sparsMethod = "threshold", 
                               thresh = 0.3,
                               dissFunc = "unsigned",
                               verbose = 3,
                               seed = 123456)

props_single <- netAnalyze(net_hembra_neg,
                           centrLCC = TRUE,
                           sPathAlgo= "automatic",
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = T, normDeg = FALSE)

summary(props_single) # summary of network properties

# Define colors
phylcol <- c("#ffa0ca", "#7570b3", "#bd0b5e", 
             "#e8b8f2", "#1b9e77", "deeppink",
             "#1f78b4", "#c51b7d", "#01665e")

# Get phyla names from the taxonomic table created before
phyla <- as.factor(gsub("p__", "", taxtab[, "Phylum"])) # sí usarlo 
names(phyla) <- taxtab[, "Phylum"] # renombrar al phylum 
phylcol_transp <- NetCoMi:::colToTransp(phylcol, 60)



fi2a <- plot(props_single, layout = "spring",
          nodeColor = "cluster",
          shortenLabels = "none",nodeTransp = 20,
          labelScale = F,
          labelLength = 15,
          nodeSize = "degree",   
          borderWidth = 1,
          borderCol = "white",
          nodeSizeSpread = 4,
          colorVec =  phylcol,
          title1 = "A) Intestinal microbial network of \nfemales fed without cry1ab/ac",      
          posCol = "gray70", labelFont = 3,
          hubLabelFont = 3, hubBorderCol = "green",
          hubBorderWidth = 5,
          negCol = "orange",
          showTitle = TRUE,
          repulsion = .84,
          cexTitle = 2, cexLabels = 1, rmSingles = "all",
          mar = c(2.5, 2.5, 4, 2.5)) 

legend(1.2, 1, cex = 1, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("gray70","orange"), 
       bty = "n", horiz = F)


### FIGURE 2B: Network female with cry1ab/ac ####
amgut_genus_pos <- phyloseq::tax_glom(physeq_hembras_cry, taxrank = "Genus")
taxtab <- amgut_genus_pos@tax_table@.Data
amgut_genus_pos@tax_table@.Data <- taxtab
rownames(amgut_genus_pos@otu_table@.Data) <- taxtab[, "Genus"]

# Network construct female whit cry1ab/ac 

net_hembra_pos <- netConstruct(amgut_genus_pos,
                               filtTax = "numbSamp",
                               filtTaxPar = list(numbSamp = 1),
                               filtSamp = "totalReads",
                               filtSampPar = list(numbTaxa = 50),
                               measure = "pearson",
                               zeroMethod = "multRepl",
                               normMethod = "clr", 
                               sparsMethod = "threshold", 
                               thresh = 0.3,
                               dissFunc = "unsigned",
                               verbose = 3,
                               seed = 123456)

props_single <- netAnalyze(net_hembra_pos,
                           centrLCC = TRUE,
                           sPathAlgo= "automatic",
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = T, normDeg = FALSE)


summary(props_single)

# Define colors
phylcol <- c("#ffa0ca", "#7570b3", "#bd0b5e", 
             "#e8b8f2", "#1b9e77", "deeppink",
             "#1f78b4", "#c51b7d", "#01665e")


p <- plot(props_single, layout = "spring",
          nodeColor = "cluster",
          shortenLabels = "none",nodeTransp = 20,
          labelScale = F,
          labelLength = 15,
          nodeSize = "degree", hubLabelFont = 3, 
          hubBorderCol = "green",  hubBorderWidth = 7, 
          borderWidth = 1,
          borderCol = "white",
          nodeSizeSpread = 4,
          colorVec =  phylcol,
          title1 = "B) Intestinal microbial network of females fed with cry1ab/ac",      
          posCol = "gray70", 
          negCol = "orange", labelFont = 3,
          showTitle = TRUE,
          repulsion = .84,
          cexTitle = 2, cexLabels =1, rmSingles = "all",
          mar = c(2.5, 2.5, 2.5, 2.5))




legend(1.2, 1, cex = 1, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("gray70","orange"), 
       bty = "n", horiz = F)


### FIGURE 3A ####
#### Compare the intestinal network between females fed with and without cry1ab / ac #####
net_season_pears <- netConstruct(data = amgut_genus_neg, 
                                 data2 = amgut_genus_pos, 
                                 filtTax = "highestVar",
                                 filtTaxPar = list(highestVar = 50),
                                 measure = "pearson", 
                                 normMethod = "clr",
                                 sparsMethod = "none", 
                                 thresh = 0.2,
                                 verbose = 3)

props_season <- netAnalyze(net_season_pears, 
                           centrLCC = FALSE, weightDeg = T,
                           avDissIgnoreInf = F,
                           sPathAlgo="automatic",
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("eigenvector"),
                           hubQuant = 0.95,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)


summary(props_season)

comp_season <- netCompare(props_season, permTest = T,
                          seed = 123, 
                          verbose = FALSE, nPerm = 100L,
                          adjust = "adaptBH")

summary(comp_season, 
        groupNames = c("cry", "Negative"),
        numbNodes = 50)

# Differential network construction
diff_season <- diffnet(net_season_pears,
                       diffMethod = "fisherTest", 
                       adjust = "lfdr")

# Differential network plot
plot(diff_season, layout = "circle",
     cexNodes = 1, labelLength = 15, labelFont = 3,
     cexLegend = .7, cexLabels = 1, labelScale = F,
     cexTitle = 1.5, borderWidth = 8,
     mar = c(2,2,3,10),edgeWidth = .5,
     edgeCol =  c("gray90", "blue", "#8462cc", "#c4ab42", "black",
                  "orange", "blue", "black", "purple"),
     legendGroupnames = c("no-cry", "cry1ab/ac"),
     title = "Differential network female",
     legendPos = c(1.3,1))



### Analysis of gut microbiota of males of Dysdercus collected on wild cotton with and without transgenes ####
### FIGURE 2C: Network male without cry1ab/ac ####
amgut_genus_neg <- phyloseq::tax_glom(physeq_macho_neg, taxrank = "Genus")
taxtab <- amgut_genus_neg@tax_table@.Data
amgut_genus_neg@tax_table@.Data <- taxtab
rownames(amgut_genus_neg@otu_table@.Data) <- taxtab[, "Genus"]

net_macho_neg <- netConstruct(amgut_genus_neg,
                              filtTax = "numbSamp",
                              filtTaxPar = list(numbSamp = 1),
                              filtSamp = "totalReads",
                              filtSampPar = list(numbTaxa = 50),
                              measure = "pearson",
                              zeroMethod = "multRepl",
                              normMethod = "clr", 
                              sparsMethod = "threshold", 
                              thresh = 0.3,
                              dissFunc = "unsigned",
                              verbose = 3,
                              seed = 123456)

props_single <- netAnalyze(net_macho_neg,
                           centrLCC = TRUE,
                           sPathAlgo= "automatic",
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = T, normDeg = FALSE)


summary(props_single)

# Define colors
phylcol <- c("#ffa0ca", "#7570b3", "#bd0b5e", 
             "#e8b8f2", "#1b9e77", "deeppink",
             "#1f78b4", "#c51b7d", "#01665e")

p <- plot(props_single, layout = "spring",
          nodeColor = "cluster",
          shortenLabels = "none",nodeTransp = 20,
          labelScale = F,
          labelLength = 15,
          nodeSize = "degree",hubLabelFont = 3,
          hubBorderCol = "green", hubBorderWidth = 7,    
          borderWidth = 1,
          borderCol = "white",
          nodeSizeSpread = 4,
          colorVec =  phylcol,
          title1 = "Intestinal microbial network of males fed without cry1ab/ac",      
          posCol = "gray70", 
          negCol = "orange",
          showTitle = TRUE, labelFont = 3,
          repulsion = .84,
          cexTitle = 2.3, cexLabels = .8, rmSingles = "all",
          mar = c(1, 3, 3, 10))

legend(1.2, 1, cex = 1, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("gray70","orange"), 
       bty = "n", horiz = F)


### FIGURE 2D: Network male with cry1ab/ac ####
amgut_genus_pos <- phyloseq::tax_glom(physeq_macho_cry, taxrank = "Genus")
taxtab <- amgut_genus_pos@tax_table@.Data
amgut_genus_pos@tax_table@.Data <- taxtab
rownames(amgut_genus_pos@otu_table@.Data) <- taxtab[, "Genus"]

net_macho_pos <- netConstruct(amgut_genus_pos,
                              filtTax = "numbSamp",
                              filtTaxPar = list(numbSamp = 1),
                              filtSamp = "totalReads",
                              filtSampPar = list(numbTaxa = 50),
                              measure = "pearson",
                              zeroMethod = "multRepl",
                              normMethod = "clr", 
                              sparsMethod = "threshold", 
                              thresh = 0.3,
                              dissFunc = "unsigned",
                              verbose = 3,
                              seed = 123456)

props_single <- netAnalyze(net_macho_pos,
                           centrLCC = TRUE,
                           sPathAlgo= "automatic",
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = T, normDeg = FALSE)

summary(props_single)

# Define colors
phylcol <- c("#ffa0ca", "deeppink", "#b399ff", 
             "#e8b8f2", "#1b9e77", "#7570b3",
             "#1f78b4", "#c51b7d", "#01665e")

p <- plot(props_single, layout = "spring",
          nodeColor = "cluster",
          shortenLabels = "none",nodeTransp = 20,
          labelScale = F, hubLabelFont = 3, 
          hubBorderCol = "green", hubBorderWidth = 7,
          labelLength = 15,
          nodeSize = "degree", 
          borderWidth = 1,
          borderCol = "white",
          nodeSizeSpread = 4,
          colorVec =  phylcol,
          title1 = "Intestinal microbial network of males fed with cry1ab/ac",      
          posCol = "gray70", labelFont = 3,
          negCol = "orange",
          showTitle = TRUE,
          repulsion = .84,
          cexTitle = 2.3, cexLabels =1, rmSingles = "all",
          mar = c(1, 3, 3, 10))

legend(1.2, 1, cex = 1, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("gray70","orange"), 
       bty = "n", horiz = F)


#### FIGURE 3B ####
#### Compare the intestinal network between males fed with and without cry1ab / ac #####
net_season_pears <- netConstruct(data = amgut_genus_neg, 
                                 data2 = amgut_genus_pos, 
                                 filtTax = "highestVar",
                                 filtTaxPar = list(highestVar = 50),
                                 measure = "pearson", 
                                 normMethod = "clr",
                                 sparsMethod = "none", 
                                 thresh = 0.2,
                                 verbose = 3)

props_season <- netAnalyze(net_season_pears, 
                           centrLCC = FALSE, weightDeg = T,
                           avDissIgnoreInf = F,
                           sPathAlgo="automatic",
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("eigenvector"),
                           hubQuant = 0.95,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)

comp_season <- netCompare(props_season, permTest = T, 
                          verbose = FALSE, nPerm = 100L,
                          adjust = "adaptBH")

summary(comp_season, 
        groupNames = c("cry", "negative"),
        numbNodes = 50, digits = 3, digitsPval = 3)


net_season_pears <- netConstruct(data = amgut_genus_neg, 
                                 data2 = amgut_genus_pos, 
                                 filtTax = "highestVar",
                                 filtTaxPar = list(highestVar = 50),
                                 measure = "pearson", 
                                 normMethod = "clr",
                                 sparsMethod = "none", 
                                 thresh = 0.2,
                                 verbose = 3)

# Differential network construction
diff_season <- diffnet(net_season_pears,
                       diffMethod = "fisherTest", 
                       adjust = "lfdr")

# Differential network plot
plot(diff_season, layout = "circle",
     cexNodes = 1, labelLength = 15, labelFont = 3,
     cexLegend = .8, cexLabels = 1, labelScale = F,
     cexTitle = 1.5, borderWidth = 8,
     mar = c(2,2,3,10),edgeWidth = .5,
     edgeCol =  c("gray90", "blue", "#8462cc", "#c4ab42", "black",
                  "orange", "blue", "black", "purple"),
     legendGroupnames = c("no-cry", "cry1ab/ac"), 
     title = "Differential network male",
     legendPos = c(1.2,1))



### Taxa dominant ####
## Phylum level ####
p0.gen <- aggregate_taxa(physeq_filtro,"Phylum")
x.d <- dominant_taxa(physeq_filtro,level = "Phylum", 
                     group="PCR_result")

x.d$dominant_overview 

ntaxa(physeq_filtro)
get_taxa_unique(physeq_filtro, "Phylum")
get_taxa_unique(physeq_hembras_cry, "Phylum")
get_taxa_unique(physeq_hembras_neg, "Phylum")

## Class level ####
p0.gen <- aggregate_taxa(physeq_filtro,"Class")
x.d <- dominant_taxa(physeq_filtro,level = "Class", 
                     group="PCR_result")

x.d$dominant_overview 

get_taxa_unique(physeq_filtro, "Class")
get_taxa_unique(physeq_hembras_cry, "Class")
get_taxa_unique(physeq_hembras_neg, "Class")

## Family level ####
p0.gen <- aggregate_taxa(physeq_filtro,"Family")
x.d <- dominant_taxa(physeq_filtro,level = "Family", 
                     group="PCR_result")

x.d$dominant_overview

get_taxa_unique(physeq_filtro, "Family")
get_taxa_unique(physeq_hembras_cry, "Family")
get_taxa_unique(physeq_hembras_neg, "Family")

##### WILCOX TEST between females fed with and without cry1ab / ac ####
(ps_clr <- microbiome::transform(physeq_hembras, 'clr'))          
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$PCR_result <- phyloseq::sample_data(ps_clr)$PCR_result

#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ PCR_result, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ PCR_result, data = df)$p.value
}

wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -PCR_result) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       

#Unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest()

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.07) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())

wilcox_results<-wilcox_results %>% filter(Phylum != "NA")
wilcox_results$p_value<- round(wilcox_results$p_value, 10)
wilcox_results$BH_FDR<- round(wilcox_results$BH_FDR, 10)

wilcox_results 

##### WILCOX TEST between males fed with and without cry1ab / ac ####

(ps_clr <- microbiome::transform(physeq_machos, "clr"))          
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$PCR_result <- phyloseq::sample_data(ps_clr)$PCR_result

#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ PCR_result, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ PCR_result, data = df)$p.value
}

wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -PCR_result) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       

#Unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest()

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.07) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())


wilcox_results<-wilcox_results %>% filter(Phylum != "NA")
wilcox_results$p_value<- round(wilcox_results$p_value, 3)
wilcox_results$BH_FDR<- round(wilcox_results$BH_FDR, 3)

wilcox_results


######## END ###########

