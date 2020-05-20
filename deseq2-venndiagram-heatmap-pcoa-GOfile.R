# # Install the bioconductor packages (you only need to do this once!)
# if (!requireNamespace("BiocManager"))
# install.packages("BiocManager")
# BiocManager::install(c("DESeq2","arrayQualityMetrics"))
# install.packages("data.table")

#BiocManager::install(c("readxl"))

library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers

# Install these from CRAN, if you haven't already
library("pheatmap") # for pretty heatmaps
library("tidyverse") # for wrangling, plotting
library("VennDiagram") # for making Venn diagram
library("vegan") # for distance matrix calculations
library("ape") # for PCoA
library("wesanderson")

pal1 <- wes_palette("Darjeeling1",1)
pal2 <- wes_palette("IsleofDogs1",1)

setwd("/Volumes/wrightlab/MCav_VibrioChallenge/RNASeq/deseq_host/") # change to your working directory

# Load data ---------------------------------------------------------------
countdata <- read.table("../allcounts_Mcav_CladeC.txt",header=TRUE,row.names = 1) 
head(countdata)
tail(countdata)
nrow(countdata)
# 36106 genes mapped

# simplify sample names
names(countdata)
names(countdata) <- gsub("MC","",names(countdata)) # replace the "MC" with nothing
names(countdata) <- gsub(".fq.sam.counts","",names(countdata)) # replace the ".fq.sam.counts" with nothing
names(countdata)

# make conditions table in excel and load it here
# or make it here if it's a simple design
conds <- read.csv("../conds.csv")
head(conds)

# VERY IMPORTANT!
# The sample names in the "conds" table must match the sample names in the counts matrix (column names)
# The order must be the same. Be sure to check it!

# First, check that the names in the countdata file exist as sample names in conds
table(names(countdata) %in% conds$sam)
# Yes! All of the names in count data are also in the conds table

# Next, check that the order matches
table(names(countdata) == conds$sam)
# No! The order is totally off. Let's fix it.

# Reorder the samples conds table to match the order in the counts matrix
conds <- conds[match(names(countdata), conds$sam),]
head(conds)

# check that it worked
table(names(countdata) == conds$sam)
# yes!

# now that everything matches, lets make our lives easier by giving our samples more informative names
head(conds)
conds$name <- paste(conds$sam,conds$bank,conds$pheno,sep="_")
head(conds)

names(countdata) <- conds$name
head(countdata)

# SEPARATE HOST AND SYMBIONT ------------------------
length(grep("Mcavernosa",row.names(countdata)))
# 21748 host genes
length(grep("isogroup",row.names(countdata)))
# 14358 symbiont genes

hostCounts <- countdata[grep("Mcavernosa",row.names(countdata)),]
head(hostCounts)
tail(hostCounts)
nrow(hostCounts)

symCounts <- countdata[grep("isogroup",row.names(countdata)),]
head(symCounts)
tail(symCounts)
nrow(symCounts)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS -------------------------------------------
# Change here depending on if you are analyzing HOST or SYMBIONT

countdata <- hostCounts
# countdata <- symCounts

#----------Total counts?
totalCounts <- colSums(countdata)
min(totalCounts) # 72909
mean(totalCounts) # 327146.5
max(totalCounts)  # 541846

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = conds,
  design = ~ pheno + bank)

save(conds, countdata, dds, file="dds_Start.Rdata")
load("dds_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(countdata,1,mean)
table(means>3)
# FALSE  TRUE 
# 13594  8154   

means3 <- names(means[means>3])
head(means3)
length(means3)
#8154

countFilt <- countdata[row.names(countdata) %in% means3,]
head(countFilt)

totalCountsFilt <- colSums(countFilt)
totalCountsFilt

min(totalCountsFilt) #71638
max(totalCountsFilt) #528713
mean(totalCountsFilt) #319547.4

# check sample order one more time... just in case!
table(names(countFilt) == as.vector(conds$name))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = countFilt,
  colData = conds,
  design = ~ pheno + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("pheno", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out2 <- c("11A_e_s", "17A_w_s", "18B_w_s", "3A_w_r", "7A_e_s")
out1 <- c("1B_w_s")

# remove out3 and out2
dim(countFilt)
countdata.out <- countFilt %>% select(-one_of(c(out2)))
dim(countdata.out)

dim(conds)
head(conds)
conds.out <- conds %>% filter(!name %in% c(out2))
dim(conds.out)

# START HERE TO ADD GROWTH DATA TO CONDS ------------
# update "conds.out" with calcification rate and tissue growth normalized to surface area
# add those two new columns to conds.out (merge function)
save(countdata.out, conds.out, file="counts_conds_filtered.RData")
load("counts_conds_filtered.RData")
load("../../tissueGrowth/growth.RData")

head(conds.out)
head(polyps)

# first, make a per-genet average of calcification rate (g/cm2/day) and normtissue (polyp/cm2/year)
avgGrowth <- polyps %>% group_by(genet) %>%
  rename(geno = genet) %>%
  summarize(meanCalc = round(mean(calcification),5),
            meanPolyp = round(mean(normtissue),2))
hist(avgGrowth$meanCalc)
hist(avgGrowth$meanPolyp)

test<- merge(conds.out, avgGrowth, by = "geno", all.x=T)
head(test)

conds.out <- test
head(conds.out)

# subset for no NA
conds.growth <- conds.out %>% filter(!is.na(meanCalc))
head(conds.growth)

head(countdata.out)
counts.growth <- countdata.out[,names(countdata.out) %in% conds.growth$name]
head(counts.growth)
dim(counts.growth)

table(names(counts.growth) %in% conds.growth$name)
table(names(counts.growth) == conds.growth$name)

test <- conds.growth[match(names(counts.growth), conds.growth$name),]
head(test)
table(names(counts.growth) == test$name)
conds.growth <- test

# look at polyp growth, try bin?
head(conds.growth)
hist(conds.growth$meanCalc,10)
test <- conds.growth %>% mutate(extremeGrowth = ifelse(meanCalc > 6e-4, "superGrowth", ifelse(meanCalc < 2e-4, "lowGrowth", "no")))
head(test)
conds.growth <- test

# Reconstruct data object (filtered and outliers removed) ------------------------------
# START HERE TO MAKE NEW MODELS THAT INCLUDE DIFFERENT CONDITIONS TO TEST -------
# change the design

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = countdata.out,
  colData = conds.out,
  design =  ~ pheno + bank)

ddsGrowth <- DESeqDataSetFromMatrix(
  countData = counts.growth,
  colData = conds.growth,
  design =  ~ extremeGrowth)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)
# Since extremeGrowth has three levels, instead of pair-wise comparisons as with pheno and bank, we'll use a LRT test
deds.growth <- DESeq(ddsGrowth, test ="LRT",reduced = ~ 1)
levels(as.factor(conds.growth$extremeGrowth))

#---Results
resPheno <- results(deds, independentFiltering = F, contrast=c("pheno","r","s"))
resPheno
# log2 fold change (MLE): pheno r vs s 
# Wald test p-value: pheno r vs s 
# here "resistant" is the numerator and "susceptible" is the denominator
# in other words, genes with a positive log2FoldChange are more highly expressed in resistant compared to susceptible genets

resBank <- results(deds, independentFiltering = F, contrast=c("bank", "e", "w"))
resBank
# log2 fold change (MLE): bank e vs w 
# Wald test p-value: bank e vs w 
# here "east" is the numerator and "west" is the denominator
# in other words, genes with a positive log2FoldChange are more highly expressed in east bank compared to west bank genets

#resPolyp <- results(deds.growth, independentFiltering = F)
#resPolyp
# log2 fold change (MLE): meanPolyp 
# Wald test p-value: meanPolyp 

resGrowth <- results(deds.growth, contrast = c ("extremeGrowth", "superGrowth", "lowGrowth"))
resGrowth
# make a new VSD file (this is a useful format for making heatmaps later)
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)
head(assay(vsd))

vsd.growth <- varianceStabilizingTransformation(ddsGrowth, blind=TRUE)
head(assay(vsd.growth))

# Write results for making heatmaps and other downstream analyses ---------------------------------------

###--------------Get pvals
head(resPheno)
valsPheno <- cbind(resPheno$pvalue, resPheno$padj)
head(valsPheno)
colnames(valsPheno)=c("pval.pheno", "padj.pheno")


head(resBank)
valsBank <- cbind(resBank$pvalue, resBank$padj)
head(valsBank)
colnames(valsBank)=c("pval.bank", "padj.bank")

#Make vsd data and pvals table
vsdpvals <- as.data.frame(cbind(assay(vsd),valsPheno, valsBank))
head(vsdpvals)
dim(vsdpvals)
# 8154   53

# for growth
head(resGrowth)
valsGrowth <- cbind(resGrowth$pvalue, resGrowth$padj)
head(valsGrowth)
colnames(valsGrowth)=c("pval.growth", "padj.growth")

vsdpvals.growth <- as.data.frame(cbind(assay(vsd.growth),valsGrowth))
head(vsdpvals.growth)
dim(vsdpvals.growth)
# 8154    43

# SAVE and LOAD (start here for downstream analyses) -------
save(avgGrowth,conds.out, countdata.out, ddsFiltOut, deds, resPheno, resBank, vsd, vsdpvals,
     conds.growth, counts.growth, ddsGrowth, deds.growth, resGrowth, vsd.growth, vsdpvals.growth,
     file="ddsCoralFiltOut.Rdata")

load("ddsCoralFiltOut.Rdata")

# Look at the results -------

# how many DE genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resPheno$padj < 0.10)
Phenop <-  resPheno [resPheno$padj <0.10,]
Phenop
#Mcavernosa08365  0.0696246192942863
#Mcavernosa23377  0.0238392770450436 ****
#Mcavernosa29127 0.00479612376999174 ****
#Mcavernosa34861  0.0835405208266223
# FALSE  TRUE 
#  8150     4

table(resBank$padj < 0.10)
Bankp <-  resBank [resBank$padj <0.10 & abs(resBank$log2FoldChange)>2,]
Bankp
# Mcavernosa04701  0.00599522500637756
# Mcavernosa23373  0.00331281889144534
# FALSE  TRUE 
#  8116    38 

#table(resPolyp$padj < 0.05)
# FALSE 
# 8154 

table(resGrowth$padj < 0.05)
# FALSE   TRUE 
# 7964    205

# how many DE genes have log2 fold changes > 2 in either direction?
table(abs(resPheno$log2FoldChange)>2)
# FALSE  TRUE 
# 8151     3 

table(abs(resBank$log2FoldChange)>2)
# FALSE  TRUE 
# 8146     8 


table(is.na(resGrowth$padj))
Growthp <-  resGrowth [resGrowth$padj <0.05 & abs(resGrowth$log2FoldChange)>4 & ! is.na(resGrowth$padj),]
Growthp
# Mcavernosa01063  0.0474005419924505
# Mcavernosa11340   0.0339241808166925
# Mcavernosa13690  0.0059095224302148
# Mcavernosa18886  0.0144139780548579
# Mcavernosa25256 0.00148694795734796
# FALSE  TRUE 
# 7707    447 

# Explore with plots ------------------------------------------------------

# Sample distance heatmap
pheatmap(cor(assay(vsd)),border_color=NA, main="SampleHeatmap")

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds, main="Dispersion Plot")
plotDispEsts(deds.growth, main="Dispersion Plot")

#MA plot
plotMA(resPheno, ylim = c(-2, 2), main="MA Plot Pheno") 
plotMA(resBank, ylim = c(-2, 2), main="MA Plot Bank")
plotMA(resGrowth, ylim = c(-2, 2), main="MA Plot Growth")

# Venn diagram ---------
sig.pheno <-  row.names(vsdpvals[vsdpvals$padj.pheno<0.05 & !is.na(vsdpvals$padj.pheno),])
sig.bank <- row.names(vsdpvals[vsdpvals$padj.bank<0.05 & !is.na(vsdpvals$padj.bank),])
sig.growth<- row.names(vsdpvals.growth[vsdpvals.growth$padj.growth<0.05 & !is.na(vsdpvals.growth$padj.growth),])

candidates <- list("Phenotype"=sig.pheno,"Bank"=sig.bank, "Growth" = sig.growth)

prettyvenn <- venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("blue", "forestgreen", "magenta"),
  alpha = 0.5,
  label.col = c("black"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("black"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)


# Make a gene expression heatmaps ------------
load("ddsCoralFiltOut.Rdata")

head(vsdpvals)
names(vsdpvals)

# subset for just the expression data
exp <- vsdpvals[c(1:49)]
head(exp)

expg <- vsdpvals.growth[c(1:41)]
head(expg)

#     Make p-value cut-offs
sig.pheno <-  row.names(vsdpvals[vsdpvals$padj.pheno<0.5 & !is.na(vsdpvals$padj.pheno),])
sig.bank <- row.names(vsdpvals[vsdpvals$padj.bank<0.05 & !is.na(vsdpvals$padj.bank),])
sig.growth <- row.names(vsdpvals[vsdpvals.growth$padj.growth<0.02 & !is.na(vsdpvals.growth$padj.growth),])

# subset expression data for significant DEGs
exp.pheno <- exp[row.names(exp) %in% sig.pheno,]
nrow(exp.pheno)

exp.bank <- exp[row.names(exp) %in% sig.bank,]
nrow(exp.bank)

exp.growth <- expg[row.names(expg) %in% sig.growth,]
nrow(exp.growth)

# Load annotations
gg <- read.delim("/Volumes/wrightlab/genomes/Mcav_genome/Mcav_gene2annot.tab", header=F)
head(gg)

# phenotype heatmap ---------
# how many are annotated?
table(row.names(exp.pheno) %in% gg$V1)

# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.pheno)){
  s=subset(gg,V1==i)
  gnames=append(gnames,paste(s$V2[1],i,sep="."))
  expg=rbind(expg,exp.pheno[i,])
} 
row.names(expg)=gnames
expl=expg
means <- apply(expl,1,mean) # means of rows
explc <- expl-means # subtracting them
head(expg)

# sort it for plotting
explc_sort <- explc %>% t() %>% as.data.frame() %>% rownames_to_column("sam") %>%
  mutate(pheno = ifelse(grepl("r", sam)==T, "resistant","susceptible"),
         bank = ifelse(grepl("w", sam)==T, "west", "east")) %>%
  arrange(pheno,bank) %>%
  column_to_rownames("sam") %>%
  select(-pheno, -bank,) %>%
  t()

# make colors 
heat.colors <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.6)(100)
annot_col <- data.frame(Phenotype = conds.out$pheno)
rownames(annot_col) <- conds.out$name
annot_color <- list(Phenotype = c(r = "coral3", s="azure3"))

# cluster plot

head(explc)
explc_sort <- explc[,c(grep("r",names(explc)),grep("s",names(explc)))]
head(explc_sort)
pheatmap(as.matrix(explc_sort), color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1], gaps_col = 13,
         clustering_distance_rows="correlation", cluster_cols=F)


pheatmap(as.matrix(explc_sort), color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1],
         clustering_distance_rows="correlation", cluster_cols=T)

`# bank heatmap ---------
means <- apply(exp.bank,1,mean) # means of rows
explc <- exp.bank-means # subtracting them

# make colors 
heat.colors <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.3)(100)

annot_col <- data.frame(Bank = conds.out$bank)
rownames(annot_col) <- names(explc)
annot_color <- list(Bank = c(e = "grey", w="black"))


# cluster plot
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1],
         clustering_distance_rows="correlation", cluster_cols=T)

# manually cluster columns based on phenotype
head(explc)
explc_sort <- explc[,c(grep("e",names(explc)),grep("w",names(explc)))]
head(explc_sort)

pheatmap(as.matrix(explc_sort),color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1],
         clustering_distance_rows="correlation", cluster_cols=F)

# growth heatmap ---------
table(row.names(exp.growth) %in% gg$V1)

# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.growth)){
  s=subset(gg,V1==i)
  gnames=append(gnames,paste(s$V2[1],i,sep="."))
  expg=rbind(expg,exp.growth[i,])
} 
row.names(expg)=gnames
expl=expg
means <- apply(expl,1,mean) # means of rows
explc <- expl-means # subtracting them
head(expg)


test <- conds.growth[match(names(expg), conds.growth$name),]
head(test)
table(names(counts.growth) == test$name)
newnames <- paste(conds.growth$name, conds.growth$extremeGrowth, sep= "_")
newnames
names(expg) = newnames

head(expg)
# sort it for plotting
explc_sort <- explc %>% t() %>% as.data.frame() %>% rownames_to_column("sam") %>%
  mutate(pheno = ifelse(grepl("r", sam)==T, "resistant","susceptible"),
         bank = ifelse(grepl("w", sam)==T, "west", "east"),
         growth = ifelse(grepl("no", sam)==T, "normal",
                         ifelse(grepl("low", sam)==T, "low growth", "super growth"))) %>%
  arrange(growth) %>%
  column_to_rownames("sam") %>%
  select(-pheno, -bank, -growth) %>%
  t()

row.names(explc_sort) <- substr(row.names(explc_sort), 1, 25)

row.names(explc_sort)
# make colors 
heat.colors <- colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.1)(100)
annot_col <- data.frame(Growth = conds.growth$extremeGrowth)
rownames(annot_col) <- names(explc)
annot_color <- list(Bank = c(e = "grey", w="black"))


# cluster plot
DEG_Growth_.1_heatmap <- pheatmap(as.matrix(explc_sort),color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1], cutree_cols = 3,
         clustering_distance_rows="correlation", cluster_cols=T)

#ggsave that plot!
ggsave(filename="DEG_Growth_.2.jpg", plot=DEG_Growth_.1_heatmap,
       width = 60, height = 40, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 200)

# manually cluster columns based on phenotype
pheatmap(as.matrix(explc_sort),color = heat.colors,cex=0.9,border_color=NA, fontsize = 16,
         annotation_col = annot_col, annotation_colors = annot_color[1],
         clustering_distance_rows="correlation", cluster_cols=F)


#Box Plots for Genes of Interest--------------- 
#Data wrangling struggles ahead :(

# BRICHOS--
#isolate the gene expression data for the genes of interest
gois <-  subset(vsdpvals[vsdpvals$padj.pheno<0.05 & !is.na(vsdpvals$padj.pheno),])
head(gois)
length(grep("23377",row.names(gois)))
length(grep("29127",row.names(gois)))

# Choose the gene you are using
row.23377 <- gois[grep("23377",row.names(gois)),]
head(row.23377)
p.23377 <- row.23377$padj.pheno[1]
row.23377 <- row.23377[c(1:49)]
colnames(row.23377)
bm.23377 <- Phenop$baseMean[2]

head(row.23377)
res_and_susc <- names(row.23377)
G23377only <- as.data.frame(cbind(t(row.23377),res_and_susc))
head(G23377only)
G23377only <- G23377only %>% mutate(pheno = ifelse(grepl("r", res_and_susc)==T, "Resistant","Susceptible"))
head(G23377only)
G23377only$Mcavernosa23377<- as.numeric(as.character(G23377only$Mcavernosa23377))

str(G23377only)

plot.23377 <- G23377only %>% ggplot(aes(x = pheno, y = Mcavernosa23377)) +
  # Make the box plot
  geom_boxplot(aes(fill=pheno)) +
  # Change the colors of the boxplot and the name in the legend
  scale_fill_manual(values=c("coral3", "azure3"),
                    name = "Resistance",
                    labels = c( "Resistant", "Susceptible")) +
  # Plot individual data points
  geom_point() +
  # Add horizontal jitter
  geom_jitter(width = 0.05) +
  # Change the colors of the points (if you want)
  scale_color_manual(values=c("black", "black")) +
  # Change the title and axis labels
  labs(title="BRICHOS Expression Levels by Phenotype", 
       x="Phenotype",
       y="Log Normalized Expression Value") + 
  # Add annotations with p-value
  annotate("text", x = 2, y = 0, size = 5,
           label = paste("p = ", format(round(p.23377,3), nsmall = 3))) +
  annotate("text", x = 1.2, y = 0, size = 5,
           label = paste("Base Mean = ", format(round(bm.23377,3), nsmall = 3))) +
  # Remove the grey background
  theme_bw() +
  # Change the sizes and placements of text and symbols
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(3,"line"))
plot.23377

#DMBT1---
gois <-  subset(vsdpvals.growth[vsdpvals.growth$padj.growth<0.05 & !is.na(vsdpvals.growth$padj.growth),])
head(gois)
length(grep("11340",row.names(gois)))


# Choose the gene you are using
row.11340 <- gois[grep("11340",row.names(gois)),]
head(row.11340)
p.11340 <- row.11340$padj.growth[1]
row.11340 <- row.11340[c(1:41)]
colnames(row.11340)
bm.11340 <- Growthp$baseMean[2]
head(row.11340)
name <- names(row.11340)
G11340only <- as.data.frame(cbind(t(row.11340),name))

head(G11340only)

data.merge <- merge(G11340only, conds.growth, by = "name")
#G11340only <- G11340only %>% mutate(pheno = ifelse(grepl("r", res_and_susc)==T, "resistant","susceptible"))
data.merge$Mcavernosa11340<- as.numeric(as.character(data.merge$Mcavernosa11340))

str(data.merge)

labs.11340 <- c("Low", "Normal", "High")

plot.11340 <- data.merge %>% ggplot(aes(x = extremeGrowth, y = Mcavernosa11340)) +
  # Make the box plot
  geom_boxplot(aes(fill=extremeGrowth)) +
  # Change the colors of the boxplot and the name in the legend
  scale_fill_manual(values=c("steelblue", "grey", "darkmagenta"),
                    name = "Growth",
                    labels = c( "Low", "Normal","High")) +
  #Label the x axis correctly
  scale_x_discrete(labels= labs.11340)+
  # Plot individual data points
  geom_point() +
  # Add horizontal jitter
  geom_jitter(width = 0.05) +
  # Change the colors of the points (if you want)
  scale_color_manual(values=c("black", "black")) +
  # Change the title and axis labels
  labs(title="DMBT1 Expression Levels by Growth Rate", 
       x="Growth Rate",
       y="Log Normalized Expression Value") + 
  # Add annotations with p-value
  annotate("text", x = 3, y = 0, size = 5,
           label = paste("p = ", format(round(p.11340,3), nsmall = 3))) +
  annotate("text", x = 1.3, y = 0, size = 5,
           label = paste("Base Mean = ", format(round(bm.11340,3), nsmall = 3))) +
  # Remove the grey background
  theme_bw() +
  # Change the sizes and placements of text and symbols
  theme(plot.title = element_text(size=20, hjust=0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(3,"line"))
plot.11340


# PCoA (for bank + pheno)---------
load("ddsCoralFiltOut.Rdata")
## I am separating bank and pheno from growth here, because a LRT test was used to generate the growth data, 
## which is not the best method for growth and pheno. 

# get variance stabilized expression data
head(vsdpvals)
names(vsdpvals)
# subset for just the expression data
exp <- vsdpvals[c(1:49)]
head(exp)

# make sure condition data match expression data
table(conds.out$name == names(exp))
# NO

# But it should be in the exp data right? 
table(conds.out$name %in% names(exp))
# Yes!

# Reorder the conds. out to match the expression data
conds.out <- conds.out[match(names(exp), conds.out$name),]
head(conds.out)

# check that it worked
table(conds.out$name == names(exp))
# yes!

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis(t(exp)~pheno+bank,
                    data=conds.out,method="manhattan")
adonisRes
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# pheno      1  10213333 10213333  1.3000 0.02657  0.051 . 
# bank       1  12777161 12777161  1.6263 0.03324  0.004 **
# Residuals 46 361407134  7856677         0.94019          
# Total     48 384397627                  1.00000     

# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors

# plotting PCoA (bank+pheno)------
margin <- .25

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA for mid by site type
plot(scores[,xaxis], scores[,yaxis],type="n", 
     main = "Coral Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.axis=1.5,
     cex.main=1.5,
     cex.lab=1.5)
# plot "spiders" connecting the samples from each bank
ordispider(scores,conds.out$pheno,label=T,cex=1)

points(scores[conds.out$pheno=="r" & conds.out$bank=="w",xaxis],
       scores[conds.out$pheno=="r" & conds.out$bank=="w",yaxis],
       col=pal1, pch=15, cex=1.5) +
  points(scores[conds.out$pheno=="r" & conds.out$bank=="e",xaxis],
         scores[conds.out$pheno=="r" & conds.out$bank=="e",yaxis],
         col=pal1, pch=16, cex=1.5) +
  points(scores[conds.out$pheno=="s" & conds.out$bank=="w",xaxis],
         scores[conds.out$pheno=="s" & conds.out$bank=="w",yaxis],
         col=pal2, pch=15, cex=1.5) +
  points(scores[conds.out$pheno=="s" & conds.out$bank=="e",xaxis],
         scores[conds.out$pheno=="s" & conds.out$bank=="e",yaxis],
         col=pal2, pch=16, cex=1.5) 

# legend of sites 
legend("topright", 
       c("Resistant", "Susceptible"),
       pch=c(8,8), 
       col=c(pal1,pal2), cex=1.5, bty = "n")
legend("bottomright", 
       c("West Bank", "East Bank"),
       pch=c(15,16), 
       col=c("grey","grey"), cex=1.5, bty = "n")

#insert p value 
legend("topleft",
       paste("Phenotype p = ",adonisRes$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=1.5, bty='n')  

legend("bottomleft", 
       paste("Bank p = ",adonisRes$aov.tab$`Pr(>F)`[2], sep=" "), 
       cex=1.5, bty='n')



# PCoA (growth)-----------------
# get variance stabilized expression data for growth
head(vsdpvals.growth)
names(vsdpvals.growth)
# subset for just the expression data
exp.g <- vsdpvals.growth[c(1:41)]
head(exp.g)
names(exp.g)
# make sure condition data match expression data
table(conds.growth$name == names(exp.g))

# compute dissimilarity indices
dd.veg.g <- vegdist(t(exp.g), "manhattan")
div.dd.veg.g <- dd.veg.g/1000
head(div.dd.veg.g)

# perform PERMANOVA  
set.seed(1)
adonisRes.g <- adonis(t(exp.g)~extremeGrowth,
                    data=conds.growth,method="manhattan")
adonisRes.g
#               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# extremeGrowth  2  23320205 11660103   1.644 0.07963  0.001 ***
# Residuals     38 269520632  7092648         0.92037           
# Total         40 292840837                  1.00000          
# 

# compute principal coordinate decomposition
dd.pcoa.g <- pcoa(div.dd.veg.g)
head(dd.pcoa.g)
scores.g <- dd.pcoa.g$vectors

# plotting PCoA (growth)----------------
plot(scores.g[,xaxis], scores.g[,yaxis],type="n", 
     main = "Coral Gene Expression",
     xlim=c(min(scores.g[,xaxis])-margin,max(scores.g[,xaxis])+margin),
     ylim=c(min(scores.g[,2])-margin,max(scores.g[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa.g$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa.g$values$Relative_eig[yaxis]*100,1),"%)",sep=""),
     cex.axis=1.5,
     cex.main=1.5,
     cex.lab=1.5)
# plot "spiders" connecting the samples from each bank
ordispider(scores.g,conds.growth$extremeGrowth,label=F, cex=.25)
# plot the supergrowth samples
points(scores.g[conds.growth$extremeGrowth=="superGrowth",xaxis],
       scores.g[conds.growth$extremeGrowth=="superGrowth",yaxis],
       col="darkmagenta", pch=15, cex=1.5) +
  # plot the lowgrowth samples
  points(scores.g[conds.growth$extremeGrowth=="lowGrowth",xaxis],
         scores.g[conds.growth$extremeGrowth=="lowGrowth",yaxis],
         col="steelblue", pch=15, cex=1.5)
  #plot the normal samples
points(scores.g[conds.growth$extremeGrowth=="no",xaxis],
       scores.g[conds.growth$extremeGrowth=="no",yaxis],
       col="grey", pch=15, cex=1.5)

# legend of sites 
legend("topright", 
       c("High Growers", "Low Growers", "Normal Range"),
       pch=c(15,15.15), 
       col=c("darkmagenta","steelblue","grey"), cex=1.5, bty = "n")

#insert p value 
legend("bottomright",
       paste("p = ",adonisRes.g$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=1.5, bty='n')  

# Write results for GO/KOG analysis by negative log pvalue-----------------------------------
load("ddsCoralFiltOut.Rdata")

# write GO results by phenotype, log-transformed signed pvalue
head(resPheno)
logs <- data.frame(cbind("gene"=row.names(resPheno),
                         "logP"=round(-log(resPheno$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resPheno$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#4159 3995 
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOPheno_logP.csv",sep=",")

# write GO results by phenotype, log-fold change
head(resPheno)
lfc <- data.frame(cbind("gene"=row.names(resPheno),
                         "LFC"=round(resPheno$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOPheno_LFC.csv",sep=",")

# write GO results by bank, log-transformed signed pvalue
head(resBank)
logs <- data.frame(cbind("gene"=row.names(resBank),
                         "logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#4116 4038 
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOBank_logP.csv",sep=",")

# write GO results by bank, log-fold change
head(resBank)
lfc <- data.frame(cbind("gene"=row.names(resBank),
                        "LFC"=round(resBank$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOBank_LFC.csv",sep=",")

# write GO results by growth rate, log-transformed signed pvalue
head(resGrowth)
logs <- data.frame(cbind("gene"=row.names(resGrowth),
                         "logP"=round(-log(resGrowth$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resGrowth$log2FoldChange<0]=-1  ##change to correct model
table(sign)
#-1     1 
#3907 4247 
logs$logP <- logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GOGrowth_logP.csv",sep=",")

# write GO results by growth rate, log-fold change
head(resGrowth)
lfc <- data.frame(cbind("gene"=row.names(resGrowth),
                        "LFC"=round(resGrowth$log2FoldChange,2)))

head(lfc)
write.table(lfc,quote=F,row.names=F,file="GOGrowth_LFC.csv",sep=",")





