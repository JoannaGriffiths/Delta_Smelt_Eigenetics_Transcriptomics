library("dplyr")
library("ggplot2")
library(tidyr)
library("bestNormalize")
library(lme4)
library(car)
library(emmeans)
library("multcomp")
library("outliers")
library("cowplot")

setwd("~/UCDavis/FCCL/Spawning_2021/CTMs")

ctm <- read.table("2021_CTM_data_trial1-120.txt", header = T)
ctm$CTM <- as.numeric(ctm$CTM)
ctm$FL <- as.numeric(ctm$FL)
ctm$DI<- as.factor(ctm$DI)
ctm$Rear_temp<- as.factor(ctm$Rear_temp)

#bestNormalize(ctm$CTM) #suggests orderNorm is best transformation
ctm$orderNorm_CTM <- predict(orderNorm(ctm$CTM))


#omit outliers:
grubbs.test(ctm$CTM, type = 10)
ctm<-subset(ctm, Fish_ID!="001H")
ctm<-subset(ctm, Fish_ID!="001R")


#######################
## ANCOVA, Fig1a
######################

LowHigh <- subset(ctm, DI=="L" | DI=="H")
ctm.model2.4 = lmer(orderNorm_CTM ~ DI*Rear_temp + FL + (1|System:DI), data = LowHigh)
Anova(ctm.model2.4, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI           28.450  1    4.01  0.005916 ** 
Rear_temp    20.498  1    5.08  0.006004 ** 
FL           66.274  1 1518.50 8.093e-16 ***
DI:Rear_temp  0.005  1    4.01  0.947002 
'''
########################



##Add in pedigree info
pedigree <- read.table("all_pedigree_meta_DI.txt", header=T)

ctm_pedigree <- merge(pedigree, ctm, by="Fish_ID", all.y = T)


##get means for each family by temperature replicate
treatment_family_CTMs <- ctm_pedigree %>%
  group_by(Rear_temp, Rep, AAFam) %>%
  summarise(mean_ctm = mean(CTM, na.rm=TRUE))

##count num of individuals per family
treatment_family_count <- ctm_pedigree %>%
  group_by(Rear_temp, Rep, AAFam) %>%
  count(AAFam)



#get means for each family by temperature
temp_family_CTMs <- ctm_pedigree %>%
  group_by(Rear_temp, AAFam) %>%
  summarise(mean_ctm = mean(CTM, na.rm=TRUE))

slopes <- temp_family_CTMs %>%
  group_by(AAFam) %>%
  filter(n_distinct(Rear_temp) > 1) %>%  # Only keep families with both temps
  summarize(slope = coef(lm(mean_ctm ~ Rear_temp))[2]) #2 keeps only slopes


family_DI <- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign/Family_DIs.txt")


slope_DI <- merge(slopes, family_DI, by="AAFam")
slope_DI$Offspring_DI <- as.numeric(slope_DI$Offspring_DI)

slope_DI_L <- subset(slope_DI, Offspring_DI < 7)
slope_DI_H <- subset(slope_DI, Offspring_DI > 9)
slope_DI_L$DI <- 'L'
slope_DI_H$DI <- 'H'
slope_DI_LH <- rbind(slope_DI_L, slope_DI_H)
slope_DI_LH$DI <- factor(slope_DI_LH$DI,levels = c("L", "H"))



##################
# CTM Slope Models, Fig1c
##################

slopes.model1 = lm(slope ~ DI, data = slope_DI_LH)
Anova(slopes.model1, test.statistic = "F")
'''
Response: slope
          Sum Sq Df F value   Pr(>F)   
DI        12.013  1  9.0801 0.003902 **
Residuals 72.763 55 
'''


######################
#Figures
######################
#https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2

## copied over from smelt_ANOVA_2021_10-3-24-LAPTOP.R
LowHigh <- subset(ctm, DI=="L" | DI=="H")
LowHigh$DI <- factor(LowHigh$DI,levels = c("L", "H"))
LowHigh$Rear_temp <- factor(LowHigh$Rear_temp,levels = c("15", "18"))

plot1 <- ggplot(data=LowHigh, aes(x=DI, y=CTM, color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))

family_DI <- read.delim2("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign/Family_DIs.txt")
temp_family_CTMs_DI <- merge(temp_family_CTMs, family_DI, by="AAFam")
temp_family_CTMs_DI$Offspring_DI <- as.numeric(temp_family_CTMs_DI$Offspring_DI)
temp_family_CTMs_DI_L <- subset(temp_family_CTMs_DI, Offspring_DI < 7)
temp_family_CTMs_DI_H <- subset(temp_family_CTMs_DI, Offspring_DI > 9)
temp_family_CTMs_DI_L$DI <- 'L'
temp_family_CTMs_DI_H$DI <- 'H'
temp_family_CTMs_DI_LH <- rbind(temp_family_CTMs_DI_L, temp_family_CTMs_DI_H)

plot2 <- ggplot(data=temp_family_CTMs_DI_LH, aes(x=Rear_temp, y=mean_ctm, group = AAFam)) +
  geom_point(size = 1, shape = 21) +
  geom_line() +
  labs(y=expression("Family Mean CTM (°C)"), x=expression("Rearing Temperature"), fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), 
        legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))


plot3 <- ggplot(data=slope_DI_LH, aes(x=DI, y=slope, color=DI)) +
  geom_boxplot() +
  geom_point(aes(fill = DI), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("grey", "black")) +
  scale_fill_manual(values=c("grey", "black")) +
  labs(y=expression("plasticity (family slope 15◦C-18◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.position="none",
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))


windows()
plot_grid(plot1, plot2, plot3, labels = "AUTO")
################




hist(treatment_family_count$n)
count_DI <- merge(treatment_family_count, family_DI, by="AAFam")

count_DI$Offspring_DI <- as.numeric(count_DI$Offspring_DI)
count_DI_L <- subset(count_DI, Offspring_DI < 7)
count_DI_H <- subset(count_DI, Offspring_DI > 9)
count_DI_L$DI <- 'L'
count_DI_H$DI <- 'H'
count_DI_LH <- rbind(count_DI_L, count_DI_H)


windows()
ggplot(data=count_DI_LH, aes(x=n, fill=DI)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c( "lightsalmon", "lightsteelblue3")) +
  labs(fill="")

windows()
ggplot(data=count_DI_LH, aes(x=n, fill=DI)) +
  geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c( "lightsalmon", "lightsteelblue3")) +
  labs(fill="")



