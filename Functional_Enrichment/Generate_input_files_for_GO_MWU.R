
'''
##Run these two lines on the bash/linux command line, not in R (removes duplicated entries)
perl nrify_GOtable.pl Htranspacficus_GOannot.txt > Htranspacficus_GOannot_nonredun.txt
sed -e "s/\r//g" Htranspacficus_GOannot_nonredun.txt > Htranspacficus_GOannot_nonredun2.txt
'''

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/RNASeq/Functional_Enrichment/GO_MWU-master/GO_MWU-master")
library("ape")
library(zoo)

##Test annotation file looks formatted correctly
goAnnotations_test <- read.delim2("Htranspacficus_GOannot_nonredun2.txt", header=F)

#######
#subset annotation file for either RNASeq or methylome background dataset
######

##rna
load("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/RNASeq/count_matrix_transformed_treatmentFilter.RData")

RNA_background <- rownames(df_log)
RNA_background <- data.frame(RNA_background)
colnames(RNA_background) <- "V1"
goAnnotations_background <- merge(RNA_background, goAnnotations_test, by="V1")
write.table(goAnnotations_background, file = "Htranspacficus_GOannot_nonredun2_RNAbackground.txt",quote=F,row.names=F, sep="\t")

##methylome temp test
meth_temp_back <- read.delim2("../../../../Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/background_no_intergenic_genenames_Temp.txt", header=T)
colnames(meth_temp_back) <- "V1"
goAnnotations_background <- merge(meth_temp_back, goAnnotations_test, by="V1")
write.table(goAnnotations_background, file = "Htranspacficus_GOannot_nonredun2_Temp_Meth_background.txt",quote=F,row.names=F, sep="\t")

##methylome DI test
meth_DI_back <- read.delim2("../../../../Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/background_no_intergenic_genenames_DI.txt", header=T)
colnames(meth_DI_back) <- "V1"
goAnnotations_background <- merge(meth_DI_back, goAnnotations_test, by="V1")
write.table(goAnnotations_background, file = "Htranspacficus_GOannot_nonredun2_DI_Meth_background.txt",quote=F,row.names=F, sep="\t")


##for rna and methylome:
colnames(goAnnotations_background) <- c("x", "GO")
goAnnotations_background$Sig <- "0"

## Read in one "input_test" file at a time and test that there is overlap in the set of significant genes/DMRs and the background set of genes (sanity check)
input_test<- read.delim2("global_Temp_sig_up.txt")
input_test<- read.delim2("global_Temp_sig_down.txt")
input_test<- read.delim2("global_DI_sig_up.txt")
input_test<- read.delim2("global_DI_sig_down.txt")
input_test<- read.delim2("global_inter_sig_up.txt")
input_test<- read.delim2("global_inter_sig_down.txt")

##Epigenetics enrichment
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/hyper_no_intergenic_genenames.txt")
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/hypo_no_intergenic_genenames.txt")

input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/hyper_no_intergenic_genenames.txt")
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/hypo_no_intergenic_genenames.txt")

input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/hyper_promoter_genenames.txt")
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/hyper_GeneBody_genenames.txt")

input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/hypo_promoter_genenames.txt")
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut/hypo_GeneBody_genenames.txt")

input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/hyper_promoter_genenames.txt")
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/hyper_GeneBody_genenames.txt")

input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/hypo_promoter_genenames.txt")
input_test<- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output/hypo_GeneBody_genenames.txt")

##for epigenetics input only:
colnames(input_test) <- c("x")


input_test$Sig <- "1"
input_test <- unique(input_test)
input_merged <- merge(goAnnotations_background, input_test, by="x", all=T)
input_merged2 <- data.frame(t(apply(input_merged, 1, na.locf, na.rm = FALSE)))

##test genes are matching
matched <- subset(input_merged2, Sig.y=="1") # yes all have a match

##test how many matched genes have GO term
matched2 <- subset(matched, Sig.x=="0") #just over half have a GO term
##533/737 for Temp hyper no intergenic
##640/932 for Temp hypo no intergenic
##409/695 for DI hyper no intergenic
##247/431 for DI hypo no intergenic
##26/51 for Temp hyper promoter
##507/686 for Temp hyper gene body
##24/68 for temp hypo promoter
##616/865 temp hypo gene body
## 31/77 for DI hyper promoter
##380/621 for DI hyper gene body
## 45/89 for DI hypo prom
## 202/342 for DI hypo gene body

##308/657 for Temp upregulated
##587/945 for Temp downregulated
##119/246 for DI upregulated
##51/107 for DI downregulated

## Write input files for GO enrichment. Select only one below that matches corresponding "input_test" file above.
write.csv(input_merged2[,c(1,4)], file="global_Temp_sig_up_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="global_Temp_sig_down_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="global_DI_sig_up_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="global_DI_sig_down_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="global_inter_sig_up_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="global_inter_sig_down_fisher.csv", quote=F,row.names=F, col.names = F)

write.csv(input_merged2[,c(1,4)], file="DMR_Temp_hyper_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="DMR_DI_hypo_fisher.csv", quote=F,row.names=F, col.names = F)

write.csv(input_merged2[,c(1,4)], file="DMR_Temp_hyper_promoter_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="DMR_Temp_hyper_GeneBody_fisher.csv", quote=F,row.names=F, col.names = F)

write.csv(input_merged2[,c(1,4)], file="DMR_Temp_hypo_promoter_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="DMR_Temp_hypo_GeneBody_fisher.csv", quote=F,row.names=F, col.names = F)

write.csv(input_merged2[,c(1,4)], file="DMR_DI_hyper_promoter_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="DMR_DI_hyper_GeneBody_fisher.csv", quote=F,row.names=F, col.names = F)

write.csv(input_merged2[,c(1,4)], file="DMR_DI_hypo_promoter_fisher.csv", quote=F,row.names=F, col.names = F)
write.csv(input_merged2[,c(1,4)], file="DMR_DI_hypo_GeneBody_fisher.csv", quote=F,row.names=F, col.names = F)

