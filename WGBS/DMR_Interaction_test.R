
library(VennDiagram)
library(dplyr)
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR/Temp_comp_outut")



##################################
##Subset dataset into low and high DI fish, then rerun DMRichR for temperature effect
##################################

library(dplyr)

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR")

##import methylation data calculated per low DI individuals for temperature effect
low_meth <- read.delim2("low_DI_temp_output/DMR_individual_smoothed_methylation.txt", header=T)
low_meth <- low_meth[,c(1,2,3,4,15)]
low_meth$region <- paste(low_meth$seqnames, low_meth$start, low_meth$end, sep = "_")
hyper_low_meth <- subset(low_meth, direction=="Hypermethylated")
hypo_low_meth <- subset(low_meth, direction=="Hypomethylated")

##import methylation data calculated per high DI individuals for temperature effect
high_meth <- read.delim2("high_DI_temp_output/DMR_individual_smoothed_methylation.txt", header=T)
high_meth <- high_meth[,c(1,2,3,4,15)]
high_meth$region <- paste(high_meth$seqnames, high_meth$start, high_meth$end, sep = "_")
hyper_high_meth <- subset(high_meth, direction=="Hypermethylated")
hypo_high_meth <- subset(high_meth, direction=="Hypomethylated")



##foverlap documentation: https://www.rdocumentation.org/packages/data.table/versions/1.15.4/topics/foverlaps
##helpful example code: https://stackoverflow.com/questions/58750440/how-to-find-overlapping-regions-between-two-data-frames-based-on-conditions
library(data.table)

hyper_low_meth <- data.table(hyper_low_meth)
hyper_high_meth <- data.table(hyper_high_meth)
hypo_low_meth <- data.table(hypo_low_meth)
hypo_high_meth <- data.table(hypo_high_meth)

##both hyper overlap
setkey(hyper_high_meth, seqnames, start, end)
hyper_overlap <- foverlaps(hyper_low_meth, hyper_high_meth, by.x=c("seqnames", "start", "end"), nomatch = 0L) #minoverlap not yet implemented, see when it will be :(
##130 regions overlap

#hypo low, hyper high
setkey(hyper_high_meth, seqnames, start, end)
hypoL_hyperH_overlap <- foverlaps(hypo_low_meth, hyper_high_meth, by.x=c("seqnames", "start", "end"), nomatch = 0L)
##0 regions overal

##hyper low, hypo high
setkey(hypo_high_meth, seqnames, start, end)
hyperL_hypoH_overlap <- foverlaps(hyper_low_meth, hypo_high_meth, by.x=c("seqnames", "start", "end"), nomatch = 0L)
##0 regions overlap

##both hypo overlap
setkey(hypo_high_meth, seqnames, start, end)
hypo_overlap <- foverlaps(hypo_low_meth, hypo_high_meth, by.x=c("seqnames", "start", "end"), nomatch = 0L) #minoverlap not yet implemented, see when it will be :(
##116 regions overlap


###############
# Venn Diagram
###############

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/Epigenetics/Analyses/DiffMethy/DMRichR")

library(VennDiagram) 

##both hyper
venn.diagram(list(Low_DI = 1:1595, High_DI = 1466:1767), fill = c("green", "yellow"), alpha = c(0.5, 0.5), lwd =0, "venn_diagram_both_hyper.tiff")

##hyper low DI, hypo high DI
venn.diagram(list(Low_DI = 1:1595, High_DI = 1596:2275), fill = c("green", "yellow"), alpha = c(0.5, 0.5), lwd =0, "venn_diagram_hyperLow_hypoHigh.tiff")

##hypo low DI, hyper high DI
venn.diagram(list(Low_DI = 1:462, High_DI = 463:764), fill = c("green", "yellow"), alpha = c(0.5, 0.5), lwd =0, "venn_diagram_hypoLow_hyperHigh.tiff")

##both hypo
venn.diagram(list(Low_DI = 1:462, High_DI = 347:1026), fill = c("green", "yellow"), alpha = c(0.5, 0.5), lwd =0, "venn_diagram_both_hypo.tiff")
