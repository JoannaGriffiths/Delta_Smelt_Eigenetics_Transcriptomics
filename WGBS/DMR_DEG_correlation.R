
#############
#Temp Effect
############
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/Temp_comp_outut")

##import methylation data calculated per individual
ind_meth <- read.delim2("DMR_individual_smoothed_methylation.txt", header=T)

#meta <- read.delim2("../../methyl_meta.txt")

##import gene names from DMRs (see R scipt for annotating DMRs output)
hyper_promoter <- read.delim2("hyper_promoter_corr_input.txt", header = T)
hyper_intron <- read.delim2("hyper_intron_corr_input.txt", header = T)
hyper_exon <- read.delim2("hyper_exon_corr_input.txt", header = T)

hyper_promoter$Region <- paste(hyper_promoter$Chr, hyper_promoter$Start, hyper_promoter$End, sep = "_")
hyper_intron$Region <- paste(hyper_intron$Chr, hyper_intron$Start, hyper_intron$End, sep = "_")
hyper_exon$Region <- paste(hyper_exon$Chr, hyper_exon$Start, hyper_exon$End, sep = "_")

hyper_all <- rbind(hyper_promoter, hyper_exon, hyper_intron)
hyper_genebody <- rbind(hyper_exon, hyper_intron)


hypo_promoter <- read.delim2("hypo_promoter_corr_input.txt", header = T)
hypo_intron <- read.delim2("hypo_intron_corr_input.txt", header = T)
hypo_exon <- read.delim2("hypo_exon_corr_input.txt", header = T)

hypo_promoter$Region <- paste(hypo_promoter$Chr, hypo_promoter$Start, hypo_promoter$End, sep = "_")
hypo_intron$Region <- paste(hypo_intron$Chr, hypo_intron$Start, hypo_intron$End, sep = "_")
hypo_exon$Region <- paste(hypo_exon$Chr, hypo_exon$Start, hypo_exon$End, sep = "_")

hypo_all <- rbind(hypo_promoter, hypo_exon, hypo_intron)
hypo_genebody <- rbind(hypo_exon, hypo_intron)


##Match with gene names with ind_meth table
ind_meth$start2 <- ind_meth$start + 1
ind_meth$Region <- paste(ind_meth$seqnames, ind_meth$start2, ind_meth$end, sep = "_")

all <- rbind(hyper_all, hypo_all)
ind_meth_gene_names <- merge(ind_meth, all, by="Region")

##change names of samples to simple number
colnames(ind_meth_gene_names) <- c("Region", "seqnames", "start", "end","width","strand","L","area","beta","stat","pval","qval","index.start", "index.end","index.width", "direction","difference",  "X307","X315",  "X317",  "X321",  "X323",  "X329",  "X331","X337",  "X339",  "X345",  "X347",  "X351",  "X353","X359",  "X361",  "367",  "X373",  "X375",  "X379","X387",  "X389",  "X417",  "X316",  "338",  "X402","X425",  "X431",  "X433",  "439",  "X451",  "457","X459",  "X465",  "X471",  "X489",  "X495",  "X508","X319",  "X325",  "X349",  "X313",  "401",  "403","409",  "415",  "302",  "308",  "310",  "318","322",  "324",  "330",  "336",  "344",  "346","350",  "358",  "360",  "366",  "372",  "374","378",  "386",  "388",  "394",  "400",  "408","410",  "414",  "416",  "418",  "427",  "477","481",  "483",  "501",  "505",  "507",  "513" ,"519",  "424",  "426",  "430",  "432",  "438","444",  "450",  "452",  "458",  "464",  "470","482",  "490",  "496",  "502",  "506",  "514","520",  "301",  "start2","Chr", "Start","End","Gene_Name")

#save(ind_meth_gene_names, file="Temp_ind_meth_gene_names.RData")
load("Temp_ind_meth_gene_names.RData")

##read in DEG for Temperature effect
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/RNASeq")
global_Temp_sig_up <- read.delim2("global_Temp_sig_up.txt", header=T)
colnames(global_Temp_sig_up) <- "Gene_Name"
global_Temp_sig_down <- read.delim2("global_Temp_sig_down.txt", header=T)
colnames(global_Temp_sig_down) <- "Gene_Name"

global_Temp_sig <- rbind(global_Temp_sig_up, global_Temp_sig_down)

##load gene count matrix that is normalized with variable name df_log
load("count_matrix_transformed_treatmentFilter.RData")
rm(y)

##filter gene count matrix with just significant genes by Temperature
library(tibble)
df_log <- data.frame(df_log)
df_log <- tibble::rownames_to_column(df_log, "Gene_Name")

Temp_df_log <- merge(global_Temp_sig, df_log, by="Gene_Name", all.x = T)
colnames(Temp_df_log) <- c("Gene_Name", "X317","X323","X347","X353","X375","X402","X417","X508","X315","X316","X331","X339","X361","X389","X433","X459","X319","X325","X349","X313","X321", "X329", "X351","X359","X379","X387","X431","X471","X495","X307", "X309","X337","X345", "X373","X425","X451","X465","X489")
save(Temp_df_log, file="Temp_df_log.RData")
load("Temp_df_log.RData")




#######################
#Correlation by gene in both DMR and DEG
######################
setwd("C:/Users/joanna.griffiths/Documents/UCDavis")
load("count_matrix_transformed_treatmentFilter.RData")
load("Temp_df_log.RData")
load("Temp_ind_meth_gene_names.RData")
load("Temp_DMR_DEG_model.RData") #output from this script

gene_rna <- as.data.frame(t(Temp_df_log))
colnames(gene_rna) <- gene_rna[1,]
gene_rna <- gene_rna[-c(1),]




##gene 1, test
gene_1_rna <- Temp_df_log[1,]
gene_1_rna <- as.data.frame(t(gene_1_rna))
gene_1_dmr <- ind_meth_gene_names[1,]
gene_1_dmr <- gene_1_dmr[,c(18:121)]
##


ind_meth_gene_names2 <- ind_meth_gene_names[,c(18:121)]
rna_dmr <- merge(Temp_df_log, ind_meth_gene_names2, by="Gene_Name")

rna_transform <- as.data.frame(t(rna_dmr[,c(1:39)]))
colnames(rna_transform) <- rna_transform[1,]
rna_transform <- rna_transform[-c(1),]
rownames(rna_transform) <- c("AWJG317", "AWJG323", "AWJG347", "AWJG353", "AWJG375", "AWJG402", "AWJG417", "AWJG508", "AWJG315", "AWJG316", "AWJG331", "AWJG339", "AWJG361", "AWJG389",
"AWJG433", "AWJG459", "AWJG319", "AWJG325", "AWJG349", "AWJG313", "AWJG321", "AWJG329", "AWJG351", "AWJG359", "AWJG379", "AWJG387", "AWJG431",
"AWJG471", "AWJG495", "AWJG307", "AWJG309", "AWJG337", "AWJG345", "AWJG373", "AWJG425", "AWJG451", "AWJG465", "AWJG489")

dmr_transform <- as.data.frame(t(rna_dmr[,c(1,40:142)]))
colnames(dmr_transform) <- dmr_transform[1,]
dmr_transform <- dmr_transform[-c(1,101:105),]

rna_dmr_transform <- merge(rna_transform, dmr_transform, by="row.names")

gene_list <- colnames(dmr_transform)


rna_transform$Indiv <- rownames(rna_transform)
dmr_transform$Indiv <- rownames(dmr_transform)

##for testing loop on one gene before testing for all genes
#gene_list = "XM_047018313.1"

x = list()
for (i in gene_list) {
  rna_transform2 <- rna_transform[,c(i,"Indiv")]
  rna_transform2$Gene <- rep(i, nrow(rna_transform2))
  
  dmr_transform2 <- dmr_transform[,c(i,"Indiv")]
  dmr_transform2$Gene <- rep(i, nrow(dmr_transform2))
  
  rna_dmr_transform2 <- merge(rna_transform2, dmr_transform2, by="Indiv")
  colnames(rna_dmr_transform2) <- c("Indiv", "RNA", "Gene.X", "DMR", "Gene.Y")
  
  x[[i]] <- rna_dmr_transform2
}

library(dplyr)
final <- bind_rows(x)
final$Gene.X <- as.factor(final$Gene.X)
final$RNA <- as.numeric(final$RNA)
final$DMR <- as.numeric(final$DMR)

#helpful link: https://forum.posit.co/t/how-to-apply-a-linear-regression-function-to-each-group-and-store-the-coefficients-in-columns/172291/2
corr_results <- final %>%
  group_by(Gene.X) %>%
  mutate(r2 = summary(lm(RNA ~ DMR))$r.squared,
         kd = summary(lm(RNA ~ DMR))$coefficients[2,1])








library(broom)
library(tidyverse)

#helpful link: https://stackoverflow.com/questions/1169539/linear-regression-and-group-by-in-r



fitted_models2 <- final %>% nest(data = -Gene.X) %>% mutate(model = map(data, ~lm(RNA~DMR, data = .)), tidied = map(model, tidy)) %>% unnest(tidied)
fitted_models2_DMR <- subset(fitted_models2, term=="DMR")
fitted_models2_DMR_sig <- subset(fitted_models2_DMR, p.value < 0.05)
fitted_models2_DMR_sig_bf <- subset(fitted_models2_DMR, p.value < (0.05/89)) #with bonferroni correction

save(final, fitted_models2, file="Temp_DMR_DEG_model.RData")
load("Temp_DMR_DEG_model.RData")

windows()
plot(rna_dmr_transform$XM_047020039.1.x, rna_dmr_transform$XM_047020039.1.y) #neg, hyper intron
windows()
plot(rna_dmr_transform$XM_047020579.1.x, rna_dmr_transform$XM_047020579.1.y) #pos, hypo exon
windows()
plot(rna_dmr_transform$XM_047020975.1.x, rna_dmr_transform$XM_047020975.1.y) #neg, hypo intron
windows()
plot(rna_dmr_transform$XM_047022316.1.x, rna_dmr_transform$XM_047022316.1.y) #pos, hyper promoter
windows()
plot(rna_dmr_transform$XM_047023175.1.x, rna_dmr_transform$XM_047023175.1.y) #neg, hypo intron
windows()
plot(rna_dmr_transform$XM_047022638.1.x, rna_dmr_transform$XM_047022638.1.y) #neg hyper intron
windows()
plot(rna_dmr_transform$XM_047027751.1.x, rna_dmr_transform$XM_047027751.1.y) #neg , hypo exon
windows()
plot(rna_dmr_transform$XM_047039291.1.x, rna_dmr_transform$XM_047039291.1.y) #neg, hypo promoter
windows()
plot(rna_dmr_transform$XM_047040555.1.x, rna_dmr_transform$XM_047040555.1.y) #neg, hyper exon
windows()
plot(rna_dmr_transform$XM_047044466.1.x, rna_dmr_transform$XM_047044466.1.y) #neg,  hypo exon
windows()
plot(rna_dmr_transform$XM_047044619.1.x, rna_dmr_transform$XM_047044619.1.y) #neg hypo intron
windows()
plot(rna_dmr_transform$XR_006957121.1.x, rna_dmr_transform$XR_006957121.1.y) #pos, hyper exon

##single test examples
rna_dmr_transform <- data.frame(sapply(rna_dmr_transform, as.numeric ))
corr_model = lm(XM_047018313.1.x ~ XM_047018313.1.y, data = rna_dmr_transform)

windows()
plot(rna_dmr_transform$XM_047018901.1.x, rna_dmr_transform$XM_047018901.1.y)


corr_model = lm(XM_047018901.1.y ~ XM_047018901.1.x, data = rna_dmr_transform)
summary <- summary(corr_model)
#https://stackoverflow.com/questions/27952653/how-to-loop-repeat-a-linear-regression-in-rel)


#############
#DI Effect
############
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_outut")

##import ind_meth data
ind_meth <- read.delim2("DMR_individual_smoothed_methylation.txt", header=T)

#meta <- read.delim2("../../methyl_meta.txt")

##import gene names from DMRs
hyper_promoter <- read.delim2("hyper_promoter_corr_input.txt", header = T)
hyper_intron <- read.delim2("hyper_intron_corr_input.txt", header = T)
hyper_exon <- read.delim2("hyper_exon_corr_input.txt", header = T)

hyper_promoter$Region <- paste(hyper_promoter$Chr, hyper_promoter$Start, hyper_promoter$End, sep = "_")
hyper_intron$Region <- paste(hyper_intron$Chr, hyper_intron$Start, hyper_intron$End, sep = "_")
hyper_exon$Region <- paste(hyper_exon$Chr, hyper_exon$Start, hyper_exon$End, sep = "_")

hyper_all <- rbind(hyper_promoter, hyper_exon, hyper_intron)
hyper_genebody <- rbind(hyper_exon, hyper_intron)


hypo_promoter <- read.delim2("hypo_promoter_corr_input.txt", header = T)
hypo_intron <- read.delim2("hypo_intron_corr_input.txt", header = T)
hypo_exon <- read.delim2("hypo_exon_corr_input.txt", header = T)

hypo_promoter$Region <- paste(hypo_promoter$Chr, hypo_promoter$Start, hypo_promoter$End, sep = "_")
hypo_intron$Region <- paste(hypo_intron$Chr, hypo_intron$Start, hypo_intron$End, sep = "_")
hypo_exon$Region <- paste(hypo_exon$Chr, hypo_exon$Start, hypo_exon$End, sep = "_")

hypo_all <- rbind(hypo_promoter, hypo_exon, hypo_intron)
hypo_genebody <- rbind(hypo_exon, hypo_intron)


##Match with gene names with ind_meth table
ind_meth$start2 <- ind_meth$start + 1
ind_meth$Region <- paste(ind_meth$seqnames, ind_meth$start2, ind_meth$end, sep = "_")

all <- rbind(hyper_all, hypo_all)
ind_meth_gene_names <- merge(ind_meth, all, by="Region")

##change names of samples to simple number
#colnames(ind_meth_gene_names) <- c("Region", "seqnames", "start", "end","width","strand","L","area","beta","stat","pval","qval","index.start", "index.end","index.width", "direction","difference",  "X307","X315",  "X317",  "X321",  "X323",  "X329",  "X331","X337",  "X339",  "X345",  "X347",  "X351",  "X353","X359",  "X361",  "367",  "X373",  "X375",  "X379","X387",  "X389",  "X417",  "X316",  "338",  "X402","X425",  "X431",  "X433",  "439",  "X451",  "457","X459",  "X465",  "X471",  "X489",  "X495",  "X508","X319",  "X325",  "X349",  "X313",  "401",  "403","409",  "415",  "302",  "308",  "310",  "318","322",  "324",  "330",  "336",  "344",  "346","350",  "358",  "360",  "366",  "372",  "374","378",  "386",  "388",  "394",  "400",  "408","410",  "414",  "416",  "418",  "427",  "477","481",  "483",  "501",  "505",  "507",  "513" ,"519",  "424",  "426",  "430",  "432",  "438","444",  "450",  "452",  "458",  "464",  "470","482",  "490",  "496",  "502",  "506",  "514","520",  "301",  "start2","Chr", "Start","End","Gene_Name")

#save(ind_meth_gene_names, file="DI_ind_meth_gene_names.RData")
load("DI_ind_meth_gene_names.RData")

##read in DEG for Temperature effect
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/RNASeq")
global_DI_sig_up <- read.delim2("global_DI_sig_up.txt", header=T)
colnames(global_DI_sig_up) <- "Gene_Name"
global_DI_sig_down <- read.delim2("global_DI_sig_down.txt", header=T)
colnames(global_DI_sig_down) <- "Gene_Name"

global_DI_sig <- rbind(global_DI_sig_up, global_DI_sig_down)

##load gene count matrix that is normalized with variable name df_log
load("count_matrix_transformed_treatmentFilter.RData")
rm(y)

##filter gene count matrix with just significant genes by Domestication
library(tibble)
df_log <- data.frame(df_log)
df_log <- tibble::rownames_to_column(df_log, "Gene_Name")

DI_df_log <- merge(global_DI_sig, df_log, by="Gene_Name", all.x = T)
colnames(DI_df_log) <- c("Gene_Name", "X317","X323","X347","X353","X375","X402","X417","X508","X315","X316","X331","X339","X361","X389","X433","X459","X319","X325","X349","X313","X321", "X329", "X351","X359","X379","X387","X431","X471","X495","X307", "X309","X337","X345", "X373","X425","X451","X465","X489")
save(DI_df_log, file="DI_df_log.RData")
load("DI_df_log.RData")




#######################
#By gene
######################
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/NOAA_computer/UCDavis/Smelt_ANOVA")
load("count_matrix_transformed_treatmentFilter.RData")
load("DI_df_log.RData")
load("DI_ind_meth_gene_names.RData")
load("DI_DMR_DEG_model.RData")

gene_rna <- as.data.frame(t(DI_df_log))
colnames(gene_rna) <- gene_rna[1,]
gene_rna <- gene_rna[-c(1),]



ind_meth_gene_names2 <- ind_meth_gene_names[,c(18:121)]
rna_dmr <- merge(DI_df_log, ind_meth_gene_names2, by="Gene_Name")

rna_transform <- as.data.frame(t(rna_dmr[,c(1:39)]))
colnames(rna_transform) <- rna_transform[1,]
rna_transform <- rna_transform[-c(1),]
rownames(rna_transform) <- c("AWJG317", "AWJG323", "AWJG347", "AWJG353", "AWJG375", "AWJG402", "AWJG417", "AWJG508", "AWJG315", "AWJG316", "AWJG331", "AWJG339", "AWJG361", "AWJG389",
                             "AWJG433", "AWJG459", "AWJG319", "AWJG325", "AWJG349", "AWJG313", "AWJG321", "AWJG329", "AWJG351", "AWJG359", "AWJG379", "AWJG387", "AWJG431",
                             "AWJG471", "AWJG495", "AWJG307", "AWJG309", "AWJG337", "AWJG345", "AWJG373", "AWJG425", "AWJG451", "AWJG465", "AWJG489")

dmr_transform <- as.data.frame(t(rna_dmr[,c(1,40:142)]))
colnames(dmr_transform) <- dmr_transform[1,]
dmr_transform <- dmr_transform[-c(1,101:105),]

rna_dmr_transform <- merge(rna_transform, dmr_transform, by="row.names")
#rna_dmr_transform <- data.frame(sapply(rna_dmr_transform, as.numeric ))

gene_list <- colnames(dmr_transform)


rna_transform$Indiv <- rownames(rna_transform)
dmr_transform$Indiv <- rownames(dmr_transform)

##for testing loop
#gene_list = "XM_047018313.1"

x = list()
for (i in gene_list) {
  rna_transform2 <- rna_transform[,c(i,"Indiv")]
  rna_transform2$Gene <- rep(i, nrow(rna_transform2))
  
  dmr_transform2 <- dmr_transform[,c(i,"Indiv")]
  dmr_transform2$Gene <- rep(i, nrow(dmr_transform2))
  
  rna_dmr_transform2 <- merge(rna_transform2, dmr_transform2, by="Indiv")
  colnames(rna_dmr_transform2) <- c("Indiv", "RNA", "Gene.X", "DMR", "Gene.Y")
  
  x[[i]] <- rna_dmr_transform2
}

library(dplyr)
#final <- bind_rows(x, .id = "column_label")
final <- bind_rows(x)
final$Gene.X <- as.factor(final$Gene.X)
final$RNA <- as.numeric(final$RNA)
final$DMR <- as.numeric(final$DMR)

#https://forum.posit.co/t/how-to-apply-a-linear-regression-function-to-each-group-and-store-the-coefficients-in-columns/172291/2
corr_results <- final %>%
  group_by(Gene.X) %>%
  mutate(r2 = summary(lm(RNA ~ DMR))$r.squared,
         kd = summary(lm(RNA ~ DMR))$coefficients[2,1])








library(broom)
library(tidyverse)

#https://stackoverflow.com/questions/1169539/linear-regression-and-group-by-in-r



fitted_models2 <- final %>% nest(data = -Gene.X) %>% mutate(model = map(data, ~lm(RNA~DMR, data = .)), tidied = map(model, tidy)) %>% unnest(tidied)
fitted_models2_DMR <- subset(fitted_models2, term=="DMR")
fitted_models2_DMR_sig <- subset(fitted_models2_DMR, p.value < 0.05) #14
fitted_models2_DMR_sig_bf <- subset(fitted_models2_DMR, p.value < (0.05/16)) # 11 with bonferroni correction

save(final, fitted_models2, file="DI_DMR_DEG_model.RData")

windows()
plot(rna_dmr_transform$XM_047020528.1.x, rna_dmr_transform$XM_047020528.1.y) #neg,  hyper intron
windows()
plot(rna_dmr_transform$XM_047025846.1.x, rna_dmr_transform$XM_047025846.1.y) #neg, hyper exon
windows()
plot(rna_dmr_transform$XM_047028340.1.x, rna_dmr_transform$XM_047028340.1.y) #neg, hyper intron
windows()
plot(rna_dmr_transform$XM_047032478.1.x, rna_dmr_transform$XM_047032478.1.y) #neg, hyper promoter
windows()
plot(rna_dmr_transform$XM_047037577.1.x, rna_dmr_transform$XM_047037577.1.y) #neg, hyper intron
windows()
plot(rna_dmr_transform$XM_047038435.1.x, rna_dmr_transform$XM_047038435.1.y) #neg , hypo intron
windows()
plot(rna_dmr_transform$XM_047038932.1.x, rna_dmr_transform$XM_047038932.1.y) #neg, hyper intron
windows()
plot(rna_dmr_transform$XM_047046573.1.x, rna_dmr_transform$XM_047046573.1.y) #neg, hyper exon
windows()
plot(rna_dmr_transform$XM_047046575.1.x, rna_dmr_transform$XM_047046575.1.y) #neg, hyper intron
windows()
plot(rna_dmr_transform$XM_047049358.1.x, rna_dmr_transform$XM_047049358.1.y) #neg,  hyper exon
windows()
plot(rna_dmr_transform$XR_006956573.1.x, rna_dmr_transform$XR_006956573.1.y) #neg, hyper promoter


##single test examples
rna_dmr_transform <- data.frame(sapply(rna_dmr_transform, as.numeric ))
corr_model = lm(XM_047018313.1.x ~ XM_047018313.1.y, data = rna_dmr_transform)

windows()
plot(rna_dmr_transform$XM_047018901.1.x, rna_dmr_transform$XM_047018901.1.y)


corr_model = lm(XM_047018901.1.y ~ XM_047018901.1.x, data = rna_dmr_transform)
summary <- summary(corr_model)
#https://stackoverflow.com/questions/27952653/how-to-loop-repeat-a-linear-regression-in-rel)



###########################
#Test for Signifcant Overlap
############################

#resource: https://www.badgrammargoodsyntax.com/compbio/2017/12/16/compbio-017-is-your-overlap-significant

rna_back <- read.delim2("Htranspacficus_GOannot_nonredun2_RNAbackground.txt")
colnames(rna_back) <- c("V1", "V2", "V3") #8425 genes
 
meth_temp_back <- read.delim2("Htranspacficus_GOannot_nonredun2_Temp_Meth_background.txt") #2571 gnes
meth_DI_back <- read.delim2("Htranspacficus_GOannot_nonredun2_DI_Meth_background.txt") #1635 genes

rna_meth_temp <- merge(rna_back, meth_temp_back, by="V1") #1302 overlap in background

rna_meth_DI <- merge(rna_back, meth_DI_back, by="V1") #659 overlap in background

##Temperature DEG_DMR overlap test
#89 DEG DMR overlap
#1244 DMRs out of 2571
#1595 DEGs out of 8425
#1302 overlap in background, total genes on platform = 1302 +(2571-1302)+(1595-1302) =2864
phyper(89-1, 1595, 8425-1302, 1244, lower.tail = FALSE, log.p = FALSE)
phyper(89-1, 1595, 2864-1595, 1244, lower.tail = FALSE, log.p = FALSE)
##pvalue=1

##DI DEG_DMR overlap test
#32 DEG DMR overlap
#965 DMRs out of 2571
#356 DEGs out of 8425
#1302 overlap in background,  total genes on platform = 1302 +(2571-1302)+(1595-1302) =2864
phyper(32-1, 356, 8425-1302, 965, lower.tail = FALSE, log.p = FALSE)
phyper(32-1, 356, 2864-1302, 965, lower.tail = FALSE, log.p = FALSE)
##pvalue=0.9925453, both not sig so don't think it matters