

# Invisible function can be used to make printed statements disappear 
invisible(lapply(c( "tidyverse", "ape", "vegan", "GGally", "adegenet", "MASS","data.table", "plyr", "lmtest", "reshape2", "Rmisc", "lmerTest","statmod"),
                 function(p){
                   if(! p %in% rownames(installed.packages())) {
                     install.packages(p)
                   }
                   library(p, character.only=TRUE)
                 }))
install.packages("rgl")

install.packages("pheatmap")

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2","edgeR","arrayQualityMetrics"))


library(DESeq2)
library(edgeR)
library(tidyverse)
library(ape)
library(vegan)
library(GGally)
library(arrayQualityMetrics)
#library(rgl)
library(adegenet)
library(MASS)
library(data.table)
library(plyr)
library(lmtest)
library(reshape2)
library(Rmisc)
library(lmerTest)
library("statmod")
library(ggplot2)
library(pheatmap)

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/RNASeq")

##Read in matrix of RSEM expected read counts 
data <- read.delim("all_quant.sf", header=TRUE, check.names= FALSE, row.names=1)

head(data)

data <- data[, c("317-R2-15-L_S19.quant", "323-R1-15-L_S11.quant", "347-R2-15-L_S20.quant", "353-R1-15-L_S3.quant", "375-R2-15-L_S32.quant", "402-R1-15-L_S12.quant", "417-R1-15-L_S40.quant", "508-R2-15-L_S33.quant", "315-R1-18-L_S8.quant", "316-R1-18-L_S16.quant", "331-R1-18-L_S6.quant", "338-R2-18-L_S5.quant", "339-R2-18-L_S22.quant", "361-R2-18-L_S23.quant", "389-R2-18-L_S24.quant", "433-R2-18-L_S39.quant", "459-R2-18-L_S38.quant", "319-R2-15-M_S34.quant", "325-R1-15-M_S13.quant", "349-R2-15-M_S35.quant", "313-R1-18-M_S30.quant", "321-R2-15-H_S7.quant", "329-R1-15-H_S9.quant", "351-R2-15-H_S31.quant", "359-R1-15-H_S25.quant", "379-R2-15-H_S17.quant", "387-R1-15-H_S26.quant",
"431-R2-15-H_S18.quant", "471-R1-15-H_S10.quant", "495-R1-15-H_S27.quant", "307-R2-18-H_S1.quant", "309-R1-18-H_S28.quant", "337-R1-18-H_S14.quant", "345-R2-18-H_S4.quant", "373-R1-18-H_S15.quant", "425-R1-18-H_S2.quant", "451-R1-18-H_S29.quant", "465-R2-18-H_S36.quant", "489-R2-18-H_S37.quant")]

##remove outlier identified on PCA
data <- subset(data, select=-`338-R2-18-L_S5.quant`)

##read in metadata
targets <- read.delim("RNA_metadata2.txt", header=TRUE, row.names=1)

##remove outlier from metadata
targets<- targets[-12,]

# Data must be rounded to nearest integer in order to be fit for negative binomial distribution
data_input <- round(data)

######################
#Plot unfiltered data
######################
# Plot distribution of unfiltered read counts across all samples, not quite unfiltered since this is already TPM
windows()
ggplot(data = data.frame(rowMeans(data_input)),
       aes(x = rowMeans.data_input.)) +
  geom_histogram(fill = "grey") +
  xlim(0, 500) +
  ylim(0, 2600) +
  theme_classic() +
  labs(title = "Distribution of unfiltered reads") +
  labs(y = "Density", x = "Raw read counts",
       title = "Read count distribution: untransformed, normalized, unfiltered")
###########################



######################
#Filtering by min expression across individuals
#########################
# Make a DGEList object for edgeR and limma voom
y <- DGEList(counts = data_input, remove.zeros = TRUE)

##Filtering by minimal expression levels per treatment group
## Keep based on presence in one treatment group
##first subset data by treatment group
y.L.15 <- y[,c(1,2,3,4,5,6,7,8)]
y.H.15 <- y[,c(17,18,19,21,22,23,24,25,26,27,28)]
y.L.18 <- y[,c(9,10,11,12,13,14,15,16)]
y.H.18 <- y[,c(20,29,30,31,32,33,34,35,36,37,38)]

colnames(y.L.15) #check correct samples are in there
y.L.15.keep <- rowSums(cpm(y.L.15)>0.5) >=8
table(y.L.15.keep)
y.L.15 <- y.L.15[y.L.15.keep,]
y.L.15.genes <- data.frame(rownames(y.L.15))
colnames(y.L.15.genes) <- "genes"

colnames(y.H.15) #check correct samples are in there
y.H.15.keep <- rowSums(cpm(y.H.15)>0.5) >=8
table(y.H.15.keep)
y.H.15 <- y.H.15[y.H.15.keep,]
y.H.15.genes <- data.frame(rownames(y.H.15))
colnames(y.H.15.genes) <- "genes"

colnames(y.L.18) #check correct samples are in there
y.L.18.keep <- rowSums(cpm(y.L.18)>0.5) >=8
table(y.L.18.keep)
y.L.18 <- y.L.18[y.L.18.keep,]
y.L.18.genes <- data.frame(rownames(y.L.18))
colnames(y.L.18.genes) <- "genes"

colnames(y.H.18) #check correct samples are in there
y.H.18.keep <- rowSums(cpm(y.H.18)>0.5) >=8
table(y.H.18.keep)
y.H.18 <- y.H.18[y.H.18.keep,]
y.H.18.genes <- data.frame(rownames(y.H.18))
colnames(y.H.18.genes) <- "genes"

all_treatment_genes <- rbind(y.L.15.genes, y.H.15.genes, y.L.18.genes, y.H.18.genes)

final_keep <- (unique(all_treatment_genes)) #left with 14,818 genes
row.names(final_keep) <- final_keep$genes
x <- merge(data_input, final_keep, by=0)
row.names(x) <- x$genes
z <- x[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39)]
y <- DGEList(counts = z, remove.zeros = TRUE)

y <- calcNormFactors(y)

# Calculate logCPM
df_log <- cpm(y, log = TRUE, prior.count = 2)

#save(data, y, df_log, file="count_matrix_transformed_treatmentFilter.RData")
#load("count_matrix_transformed_treatmentFilter.RData")
write.table(df_log, file="count_matrix_dflog.txt", quote=F, sep = "\t")
##############################################




#####################
#Plotting filtered data
######################
# Plot distribution of filtered logCPM values
windows()
ggplot(data = data.frame(rowMeans(df_log)), 
       aes(x = rowMeans.df_log.) ) +
  geom_histogram(fill = "grey") +
  theme_classic() +
  labs(y = "Density", x = "Filtered read counts (logCPM)",
       title = "Distribution of normalized, filtered read counts")
###################################


######################
#PCOA of filtered data
#####################
# Export pcoa loadings
dds.pcoa = pcoa(vegdist(t(df_log <- cpm(y, log = TRUE, prior.count = 2)),
                        method = "euclidean") / 1000)

# Create df of MDS vector loading
scores <- dds.pcoa$vectors

## Plot pcoa loadings of each sample, grouped by Temp and DI

# Calculate % variation explained by each eigenvector
percent <- dds.pcoa$values$Eigenvalues
cumulative_percent_variance <- (percent / sum( percent)) * 100


ordispider(scores, as.factor(targets$Replicates1), label = F) 
# Vectors connecting samples in same pCO2 x time group

ordilabel(scores, cex = 0.5) # Label sample IDs

# Prepare information for pcoa plot, then plot
# For Temp
color <- c("steelblue1", "tomato1")
color <- c("lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightsalmon","lightsalmon", "lightsalmon","lightsalmon","lightsalmon","lightsalmon","lightsalmon","lightsalmon", "steelblue3","steelblue3","steelblue3","salmon3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3",  "salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3")

par(mfrow = c(1, 1))
#par(xpd=TRUE)
#par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(
  scores[, 1],
  scores[, 2],
  cex = 1,
  cex.axis = 1,
  cex.lab = 1.25,
  col=color,
  pch = 19,
  xlab = paste("PC1, ", round(cumulative_percent_variance[1], 2), "%"),
  ylab = paste("PC2, ", round(cumulative_percent_variance[2], 2), "%")
)
legend("topright", legend=c("Low DI 15C", "High DI 15C", "Low DI 18C", "High DI 18C"), pch=20,col=c("lightblue2", "steelblue3", "lightsalmon", "salmon3"), cex=1)



## Add visual groupings to pcoa plot
#ordispider(scores,as.factor(targets$Temp), col=color)
targets$grouping <- factor(paste(targets$Temp, targets$Simple_DI, sep = "_"))

ordiellipse(
  scores,
  as.factor(targets$grouping),
  border = T,
  lty = 1,
  lwd = 2,
  label = F,
  col = c("lightblue2","steelblue3","lightsalmon","salmon3"),
  draw = "lines",
  alpha = 100,
  cex = 1
)




##############
#sum to zero parametrization (following limma voom manual)

#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.montana.edu/rotella/documents/502/DesignMatricesR.pdf
#https://faculty.nps.edu/sebuttre/home/r/contrasts.html
#examples: https://github.com/mariestrader/S.purp_RRBS_RNAseq_2019/blob/master/scripts/edgeR_RRBS_final.R
##################

Temp<- factor(targets$Temp, levels=c("15","18"))
Num_DI <- factor(targets$Simple_DI, levels=c("1","3"))
System <- as.numeric(targets$System)
Replicate <- as.numeric(targets$Replicate)

contrasts(Num_DI) <- contr.sum(2)
contrasts(Temp) <- contr.sum(2)
#design <- model.matrix(~Num_DI*Temp)
design <- model.matrix(~Num_DI*Temp + (1|System)) #final design for following analyses
#design <- model.matrix(~Num_DI*Temp + (1|Replicate))

Voom <- voom(y, design, plot = F)
# Fit using Voom
lm_Voom_fit <- lmFit(Voom, design)

#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models. This will help with coefficient contrast set up for topTable:
head(coef(lm_Voom_fit))

####Single contrast tests (doesn't correct for multiple testing)
#DI main effect
cont_DI <- contrasts.fit(lm_Voom_fit, coef = 2)
# Temperature main effect
cont_Temp <- contrasts.fit(lm_Voom_fit, coef = 3)
# interaction of DI and temperature
cont_inter <- contrasts.fit(lm_Voom_fit, coef=5)

cont_DI <- eBayes(cont_DI)
cont_Temp <- eBayes(cont_Temp)
cont_inter <- eBayes(cont_inter)

fit2 <- contrasts.fit(lm_Voom_fit, design)
fit2 <- eBayes(fit2)
results_mult2 <- decideTests(fit2, method="global", adjust.method = "fdr")
summary(results_mult)

DI_results <-topTable(cont_DI,
                            #coef = 2,
                            adjust.method = "fdr",
                            n = Inf)
Temp_results <-topTable(cont_Temp,
                           # coef = 3,
                            adjust.method = "fdr",
                            n = Inf)
inter_results <-topTable(cont_inter,
                        # coef = 4,
                         adjust.method = "fdr",
                         n = Inf)

# How many DEGs are associated with Temp
length(which(DI_results$adj.P.Val < 0.05))
# How many DEGs are associated with Temp
length(which(Temp_results$adj.P.Val < 0.05))
# How many DEGs are associated with Temp:Num_DI?
length(which(inter_results$adj.P.Val < 0.05))

##Write out csv files that include a table of genes with corresponding p-values
write.csv(DI_results, file = "sumzero_dif.gene.DImainEffect_treatmentFilter.csv")
write.csv(Temp_results, file = "sumzero_dif.gene.TempMainEffect_treatmentFilter.csv")
write.csv(inter_results, file = "sumzero_dif.gene.InterMainEffect_treatmentFilter.csv")

##### Same design matrix as above but now with multiple testing. These results presented in manuscript
##Make contrasts with matrix
cont.matrix <- cbind(DIMainEffect=c(0,1,0,0,0), TempMainEffect=c(0,0,1,0,0), Inter=c(0,0,0,0,1))


fit2 <- contrasts.fit(lm_Voom_fit, cont.matrix)
fit2 <- eBayes(fit2)

results_mult <- decideTests(fit2, method="global", adjust.method = "fdr")
summary(results_mult)
#plotMD(results_mult)
#vennDiagram(results_mult)
heatDiagram(results_mult, fit2$coefficients)
heatDiagram(results_mult, fit2$coefficients, primary = "TempMainEffect")

write.fit(fit2, adjust="fdr", method="global", file="globalDEresults_treatmentFilter.txt")


globalresults<-read.delim2(file="globalDEresults_treatmentFilter.txt", header=T)
##DI main effect
length(which(globalresults$P.value.adj.DIMainEffect < 0.05)) #353
##Temp main effect
length(which(globalresults$P.value.adj.TempMainEffect < 0.05)) #1602
##DI main effect
length(which(globalresults$P.value.adj.Inter < 0.05)) #27


DI_results <-topTable(fit2,
                      #coef = 2,
                      adjust.method = "fdr",
                      n = Inf)
Temp_results <-topTable(cont_Temp,
                        # coef = 3,
                        adjust.method = "fdr",
                        n = Inf)
inter_results <-topTable(cont_inter,
                         # coef = 4,
                         adjust.method = "fdr",
                         n = Inf)

####Make input for MWU Misha Matz functional enrichment
globalresults <- read.delim2("globalDEresults_treatmentFilter.txt")

global_DI_sig_up <- subset(globalresults, P.value.adj.DIMainEffect < 0.05 & t.DIMainEffect < 0)[,1]
write.table(global_DI_sig_up, file="global_DI_sig_up.txt",quote=F,row.names=F) 

global_DI_sig_down <- subset(globalresults, P.value.adj.DIMainEffect < 0.05 & t.DIMainEffect > 0)[,1]
write.table(global_DI_sig_down, file="global_DI_sig_down.txt",quote=F,row.names=F)

global_Temp_sig_up <- subset(globalresults, P.value.adj.TempMainEffect < 0.05 & t.TempMainEffect < 0)[,1]
write.table(global_Temp_sig_up, file="global_Temp_sig_up.txt",quote=F,row.names=F)

global_Temp_sig_down <- subset(globalresults, P.value.adj.TempMainEffect < 0.05 & t.TempMainEffect > 0)[,1]
write.table(global_Temp_sig_down, file="global_Temp_sig_down.txt",quote=F,row.names=F)

global_Inter_sig_up <- subset(globalresults, P.value.adj.Inter < 0.05 & t.Inter < 0)[,1]
write.table(global_Inter_sig_up, file="global_Inter_sig_up.txt",quote=F,row.names=F)

global_Inter_sig_down <- subset(globalresults, P.value.adj.Inter < 0.05 & t.Inter > 0)[,1]
write.table(global_Inter_sig_down, file="global_Inter_sig_down.txt",quote=F,row.names=F)


###################################
## heatmap 
##################################

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/RNASeq")
library(gplots) 
#install.packages("RColorBrewer")
library("RColorBrewer") 
library( "genefilter" )

load(file="count_matrix_transformed_treatmentFilter.RData")
globalresults <- read.delim2("globalDEresults_treatmentFilter.txt")

global_DI_sig <- subset(globalresults, P.value.adj.DIMainEffect < 0.05)[,1]
global_Temp_sig <- subset(globalresults, P.value.adj.TempMainEffect < 0.05)[,1]
global_inter_sig <- subset(globalresults, P.value.adj.Inter < 0.05)[,1]

##check which is up and down regulated in reference to. The results of this show that these 920 genes are upregulated at 15C. 103 gense upregulated for Low DI.
#global_Temp_sig_up <- read.delim2("global_Temp_sig_up.txt")
#global_Temp_sig <- global_Temp_sig_up[,1]
#global_Temp_sig_down <- read.delim2("global_Temp_sig_down.txt")
#global_Temp_sig <- global_Temp_sig_down[,1]
#global_DI_sig_up <- read.delim2("global_DI_sig_up.txt")
#global_DI_sig <- global_DI_sig_up[,1]

df_log_DIorder <- df_log[,c("323-R1-15-L_S11.quant",	"353-R1-15-L_S3.quant",	"402-R1-15-L_S12.quant",	"417-R1-15-L_S40.quant",	"317-R2-15-L_S19.quant",	"508-R2-15-L_S33.quant",	"375-R2-15-L_S32.quant",	"347-R2-15-L_S20.quant",	"316-R1-18-L_S16.quant",	"315-R1-18-L_S8.quant",	"331-R1-18-L_S6.quant",	"361-R2-18-L_S23.quant",	"459-R2-18-L_S38.quant",	"389-R2-18-L_S24.quant",	"433-R2-18-L_S39.quant",	"339-R2-18-L_S22.quant",	"325-R1-15-M_S13.quant",	"387-R1-15-H_S26.quant",	"495-R1-15-H_S27.quant",	"359-R1-15-H_S25.quant",	"471-R1-15-H_S10.quant",	"329-R1-15-H_S9.quant",	"349-R2-15-M_S35.quant",	"319-R2-15-M_S34.quant",	"351-R2-15-H_S31.quant",	"431-R2-15-H_S18.quant",	"321-R2-15-H_S7.quant",	"379-R2-15-H_S17.quant",	"313-R1-18-M_S30.quant",	"373-R1-18-H_S15.quant",	"425-R1-18-H_S2.quant",	"451-R1-18-H_S29.quant",	"337-R1-18-H_S14.quant",	"309-R1-18-H_S28.quant",	"465-R2-18-H_S36.quant",	"489-R2-18-H_S37.quant",	"307-R2-18-H_S1.quant",	"345-R2-18-H_S4.quant")]
col_colors<- c("lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightsalmon","lightsalmon", "lightsalmon","lightsalmon","lightsalmon","lightsalmon","lightsalmon","lightsalmon", "steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3",  "salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3")

df_log_DIorder.sorted <- sort(df_log_DIorder, decreasing=T)


par(mar=c(7,4,4,2)+1) 
png(filename='global_diff.gene.DI.heatmap2.png', width=800, height=750)

heatmap.2((df_log_DIorder)[global_DI_sig, ], scale="row", Colv = FALSE, 
          trace="none",  dendrogram="none", key=TRUE, keysize = 0.65, margins =c(10,11),
          density.info = "none",
          ColSideColors=col_colors,
          #col = colorpanel(75,"dodgerblue2","black","yellow"))
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"))

par(lend = 1)           # square line ends for the color legend
graphics.off()

##beuatiful re-ordering thanks to this site: https://www.biostars.org/p/274509/
## uses use 1-Pearson's correlation distance


df_log_TempOrder <- df_log[,c("323-R1-15-L_S11.quant", "353-R1-15-L_S3.quant", "402-R1-15-L_S12.quant",	"417-R1-15-L_S40.quant", "317-R2-15-L_S19.quant",	"508-R2-15-L_S33.quant",	"375-R2-15-L_S32.quant",	"347-R2-15-L_S20.quant",	"325-R1-15-M_S13.quant",	"387-R1-15-H_S26.quant",	"495-R1-15-H_S27.quant",	"359-R1-15-H_S25.quant", "471-R1-15-H_S10.quant",	"329-R1-15-H_S9.quant",	"349-R2-15-M_S35.quant",	"319-R2-15-M_S34.quant",	"351-R2-15-H_S31.quant",	"431-R2-15-H_S18.quant",	"321-R2-15-H_S7.quant",	"379-R2-15-H_S17.quant",	"316-R1-18-L_S16.quant", "315-R1-18-L_S8.quant",	"331-R1-18-L_S6.quant",	"361-R2-18-L_S23.quant",	"459-R2-18-L_S38.quant",	"389-R2-18-L_S24.quant",	"433-R2-18-L_S39.quant", "339-R2-18-L_S22.quant",	"313-R1-18-M_S30.quant",	"373-R1-18-H_S15.quant",	"425-R1-18-H_S2.quant",	"451-R1-18-H_S29.quant",	"337-R1-18-H_S14.quant",	"309-R1-18-H_S28.quant",	"465-R2-18-H_S36.quant",	"489-R2-18-H_S37.quant",	"307-R2-18-H_S1.quant",	"345-R2-18-H_S4.quant")]
col_colors<- c("lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2","lightblue2", "steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3", "lightsalmon","lightsalmon", "lightsalmon","lightsalmon","lightsalmon","lightsalmon","lightsalmon","lightsalmon", "salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3","salmon3")
                         
par(mar=c(7,4,4,2)+1) 
png(filename='global_diff.gene.Temp.heatmap2.png', width=800, height=750)

heatmap.2((df_log_TempOrder)[global_Temp_sig, ], scale="row", Colv = FALSE,
          trace="none", dendrogram="none", key=TRUE, keysize = 0.65, #key size controls size of legend
          margins =c(10,11),
          density.info = "none",
          ColSideColors=col_colors,
          #col = colorpanel(75,"dodgerblue2","black","yellow"))
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"))
par(lend = 1)           # square line ends for the color legend

graphics.off()




par(mar=c(7,4,4,2)+1) 
png(filename='global_diff.gene.inter.heatmap2.png', width=800, height=750)
heatmap.2((df_log_TempOrder)[global_inter_sig, ], scale="row", Colv = FALSE,
          trace="none", dendrogram="none", key=TRUE, keysize = 0.65, margins =c(10,11), density.info = "none", ColSideColors=col_colors,
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"))
par(lend = 1)           # square line ends for the color legend

graphics.off()


par(mar=c(7,4,4,2)+1) 
png(filename='test.inter.heatmap2.png', width=800, height=750)
heatmap.2((df_log_TempOrder)[global_inter_sig, ], scale="column", Colv = FALSE,
          trace="none", dendrogram="none", key=TRUE, keysize = 0.65, margins =c(10,11), density.info = "none", ColSideColors=col_colors,
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x, method="complete"))
par(lend = 1)           # square line ends for the color legend

graphics.off()




