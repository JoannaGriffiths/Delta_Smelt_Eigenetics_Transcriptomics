
library("stringr")
library(tidyr)

###########
#Temp Effect output from HOMER
##########
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut")

hyper <- read.delim2("HOMER_annotation_hyper", header = T)

hyper_promoter <- hyper[str_detect(hyper$Annotation, "promoter"), ]
hyper_intron <- hyper[str_detect(hyper$Annotation, "intron"), ]
hyper_exon <- hyper[str_detect(hyper$Annotation, "exon"), ]


hyper_promoter <- hyper_promoter %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hyper_promoter$Gene_Name <- gsub(")", "", hyper_promoter$Gene_Name)

hyper_intron <- hyper_intron %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hyper_intron <- hyper_intron %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hyper_exon <- hyper_exon %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hyper_exon <- hyper_exon %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hyper_GeneBody <- rbind(hyper_exon, hyper_intron)


write.table(hyper_promoter[,9], file="hyper_promoter_genenames.txt", sep = "\t", row.names = F)
write.table(hyper_intron[,9], file="hyper_intron_genenames.txt", sep = "\t", row.names = F)
write.table(hyper_exon[,9], file="hyper_exon_genenames.txt", sep = "\t", row.names = F)
write.table(hyper_GeneBody[,9], file="hyper_GeneBody_genenames.txt", sep = "\t", row.names = F)

hyper_promoter <- hyper_promoter[,9]
hyper_GeneBody <- hyper_GeneBody[,9]
hyper_no_intergenic <- rbind(hyper_promoter, hyper_GeneBody)

write.table(hyper_no_intergenic, file="hyper_no_intergenic_genenames.txt",sep = "\t", row.names = F)


##for DMR-RNA correlation analysis
write.table(hyper_promoter[,c(2,3,4,9)], file="hyper_promoter_corr_input.txt", sep = "\t", row.names = F)
write.table(hyper_exon[,c(2,3,4,9)], file="hyper_exon_corr_input.txt", sep = "\t", row.names = F)
write.table(hyper_intron[,c(2,3,4,9)], file="hyper_intron_corr_input.txt", sep = "\t", row.names = F)


hypo <- read.delim2("HOMER_annotation_hypo", header = T)

hypo_promoter <- hypo[str_detect(hypo$Annotation, "promoter"), ]
hypo_intron <- hypo[str_detect(hypo$Annotation, "intron"), ]
hypo_exon <- hypo[str_detect(hypo$Annotation, "exon"), ]


hypo_promoter <- hypo_promoter %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hypo_promoter$Gene_Name <- gsub(")", "", hypo_promoter$Gene_Name)

hypo_intron <- hypo_intron %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hypo_intron <- hypo_intron %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hypo_exon <- hypo_exon %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hypo_exon <- hypo_exon %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hypo_GeneBody <- rbind(hypo_exon, hypo_intron)


write.table(hypo_promoter[,9], file="hypo_promoter_genenames.txt", sep = "\t", row.names = F)
write.table(hypo_intron[,9], file="hypo_intron_genenames.txt", sep = "\t", row.names = F)
write.table(hypo_exon[,9], file="hypo_exon_genenames.txt", sep = "\t", row.names = F)
write.table(hypo_GeneBody[,9], file="hypo_GeneBody_genenames.txt", sep = "\t", row.names = F)

hypo_promoter <- hypo_promoter[,9]
hypo_GeneBody <- hypo_GeneBody[,9]
hypo_no_intergenic <- rbind(hypo_promoter, hypo_GeneBody)

write.table(hypo_no_intergenic, file="hypo_no_intergenic_genenames.txt",sep = "\t", row.names = F)


##for DMR-RNA correlation analysis
write.table(hypo_promoter[,c(2,3,4,9)], file="hypo_promoter_corr_input.txt", sep = "\t", row.names = F)
write.table(hypo_exon[,c(2,3,4,9)], file="hypo_exon_corr_input.txt", sep = "\t", row.names = F)
write.table(hypo_intron[,c(2,3,4,9)], file="hypo_intron_corr_input.txt", sep = "\t", row.names = F)

###########
#DI Effect output from HOMER
##########
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output")

hyper <- read.delim2("HOMER_annotation_hyper", header = T)

hyper_promoter <- hyper[str_detect(hyper$Annotation, "promoter"), ]
hyper_intron <- hyper[str_detect(hyper$Annotation, "intron"), ]
hyper_exon <- hyper[str_detect(hyper$Annotation, "exon"), ]


hyper_promoter <- hyper_promoter %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hyper_promoter$Gene_Name <- gsub(")", "", hyper_promoter$Gene_Name)

hyper_intron <- hyper_intron %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hyper_intron <- hyper_intron %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hyper_exon <- hyper_exon %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hyper_exon <- hyper_exon %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hyper_GeneBody <- rbind(hyper_exon, hyper_intron)


write.table(hyper_promoter[,9], file="hyper_promoter_genenames.txt", sep = "\t", row.names = F)
write.table(hyper_intron[,9], file="hyper_intron_genenames.txt", sep = "\t", row.names = F)
write.table(hyper_exon[,9], file="hyper_exon_genenames.txt", sep = "\t", row.names = F)
write.table(hyper_GeneBody[,9], file="hyper_GeneBody_genenames.txt", sep = "\t", row.names = F)

hyper_promoter <- hyper_promoter[,9]
hyper_GeneBody <- hyper_GeneBody[,9]
hyper_no_intergenic <- rbind(hyper_promoter, hyper_GeneBody)

write.table(hyper_no_intergenic, file="hyper_no_intergenic_genenames.txt",sep = "\t", row.names = F)

##for DMR-RNA correlation analysis
write.table(hyper_promoter[,c(2,3,4,9)], file="hyper_promoter_corr_input.txt", sep = "\t", row.names = F)
write.table(hyper_exon[,c(2,3,4,9)], file="hyper_exon_corr_input.txt", sep = "\t", row.names = F)
write.table(hyper_intron[,c(2,3,4,9)], file="hyper_intron_corr_input.txt", sep = "\t", row.names = F)

hypo <- read.delim2("HOMER_annotation_hypo", header = T)

hypo_promoter <- hypo[str_detect(hypo$Annotation, "promoter"), ]
hypo_intron <- hypo[str_detect(hypo$Annotation, "intron"), ]
hypo_exon <- hypo[str_detect(hypo$Annotation, "exon"), ]
hypo_GeneBody <- rbind(hyper_exon, hyper_intron)

hypo_promoter <- hypo_promoter %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hypo_promoter$Gene_Name <- gsub(")", "", hypo_promoter$Gene_Name)

hypo_intron <- hypo_intron %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hypo_intron <- hypo_intron %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hypo_exon <- hypo_exon %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
hypo_exon <- hypo_exon %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

hypo_GeneBody <- rbind(hypo_exon, hypo_intron)


write.table(hypo_promoter[,9], file="hypo_promoter_genenames.txt", sep = "\t", row.names = F)
write.table(hypo_intron[,9], file="hypo_intron_genenames.txt", sep = "\t", row.names = F)
write.table(hypo_exon[,9], file="hypo_exon_genenames.txt", sep = "\t", row.names = F)
write.table(hypo_GeneBody[,9], file="hypo_GeneBody_genenames.txt", sep = "\t", row.names = F)

hypo_promoter <- hypo_promoter[,9]
hypo_GeneBody <- hypo_GeneBody[,9]
hypo_no_intergenic <- rbind(hypo_promoter, hypo_GeneBody)

write.table(hypo_no_intergenic, file="hypo_no_intergenic_genenames.txt",sep = "\t", row.names = F)

##for DMR-RNA correlation analysis
write.table(hypo_promoter[,c(2,3,4,9)], file="hypo_promoter_corr_input.txt", sep = "\t", row.names = F)
write.table(hypo_exon[,c(2,3,4,9)], file="hypo_exon_corr_input.txt", sep = "\t", row.names = F)
write.table(hypo_intron[,c(2,3,4,9)], file="hypo_intron_corr_input.txt", sep = "\t", row.names = F)


###########
#Background set from Temperature Effect
##########
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/temp_comp_outut")

back <- read.delim2("background_HOMER_annotation", header = T)

back_promoter <- back[str_detect(back$Annotation, "promoter"), ]
back_intron <- back[str_detect(back$Annotation, "intron"), ]
back_exon <- back[str_detect(back$Annotation, "exon"), ]


back_promoter <- back_promoter %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
back_promoter$Gene_Name <- gsub(")", "", back_promoter$Gene_Name)

back_intron <- back_intron %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
back_intron <- back_intron %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

back_exon <- back_exon %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
back_exon <- back_exon %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

back_GeneBody <- rbind(back_exon, back_intron)


back_promoter <- back_promoter[,9]
back_GeneBody <- back_GeneBody[,9]
back_no_intergenic <- rbind(back_promoter, back_GeneBody)

write.table(back_no_intergenic, file="background_no_intergenic_genenames_Temp.txt",sep = "\t", row.names = F)


###########
#Background set from DI Effect
##########
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/DI_comp_output")

back <- read.delim2("background_HOMER_annotation", header = T)

back_promoter <- back[str_detect(back$Annotation, "promoter"), ]
back_intron <- back[str_detect(back$Annotation, "intron"), ]
back_exon <- back[str_detect(back$Annotation, "exon"), ]


back_promoter <- back_promoter %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
back_promoter$Gene_Name <- gsub(")", "", back_promoter$Gene_Name)

back_intron <- back_intron %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
back_intron <- back_intron %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

back_exon <- back_exon %>% separate_wider_delim(Annotation, "(", names =c("Annotation", "Gene_Name"))
back_exon <- back_exon %>% separate_wider_delim(Gene_Name, ",", names =c("Gene_Name", "number"))

back_GeneBody <- rbind(back_exon, back_intron)


back_promoter <- back_promoter[,9]
back_GeneBody <- back_GeneBody[,9]
back_no_intergenic <- rbind(back_promoter, back_GeneBody)

write.table(back_no_intergenic, file="background_no_intergenic_genenames_DI.txt",sep = "\t", row.names = F)

