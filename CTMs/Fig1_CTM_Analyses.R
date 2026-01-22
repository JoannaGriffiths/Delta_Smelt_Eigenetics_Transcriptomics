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
## ANCOVA, Fig2a
######################

LowHigh <- subset(ctm, DI=="L" | DI=="H")

##Add in pedigree info
pedigree <- read.table("all_pedigree_meta_DI.txt", header=T)
family_DI <- read.delim2("Family_DIs.txt")
ctm_pedigree <- merge(pedigree, LowHigh, by="Fish_ID", all.y = T)
ctm_pedigree_DI <- merge(ctm_pedigree, family_DI, by="AAFam")

LowHigh_ctm_pedigree_DI <- subset(ctm_pedigree_DI, assigned_DI=="L" | assigned_DI=="H")
ctm.model2.4 = lmer(orderNorm_CTM ~ assigned_DI*Rear_temp + FL + (1|System:DI) + (1|Trial), data = LowHigh_ctm_pedigree_DI)
Anova(ctm.model2.4, test.statistic = "F")
'''
                           F Df  Df.res    Pr(>F)    
assigned_DI           23.6572  1    4.84  0.005022 ** 
Rear_temp              9.9948  1    5.53  0.021816 *  
FL                    61.9056  1 1156.09 8.232e-15 ***
assign
'''

#######################
## Fig2b
######################

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
slope_DI_H <- subset(slope_DI, Offspring_DI > 10)
slope_DI_L$DI <- 'L'
slope_DI_H$DI <- 'H'
slope_DI_LH <- rbind(slope_DI_L, slope_DI_H)
slope_DI_LH$DI <- factor(slope_DI_LH$DI,levels = c("L", "H"))



##################
# CTM Slope Models, Fig2c
##################

slopes.model1 = lm(slope ~ DI, data = slope_DI_LH)
Anova(slopes.model1, test.statistic = "F")
''' ## with family greater than 0, 18 low families, and 27 high families
Response: slope
          Sum Sq Df F value  Pr(>F)  
DI         6.832  1  3.8334 0.05675 .
Residuals 76.631 43    
'''

''' ## with family greater than 1, 12 families for low, 20 families for high
           Sum Sq Df F value Pr(>F)
DI         2.4059  1  2.7498 0.1077
Residuals 26.2486 30  
'''


######################
#Plots
######################
#https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2


#########
#Figure S1
###################### numerical DI per categorical DI
ctm_pedigree_DI <- merge(family_DI, ctm_pedigree, by="AAFam")
ctm_pedigree_DI$Offspring_DI.x <- as.numeric(ctm_pedigree_DI$Offspring_DI.x)
ctm_pedigree_DI$assigned_DI <- factor(ctm_pedigree_DI$assigned_DI,levels = c("L", "M", "H"))

windows()
ggplot(data=ctm_pedigree_DI, aes(x=assigned_DI, y=Offspring_DI.x)) +
  geom_boxplot() +
  geom_point(aes(fill = assigned_DI), size = 1, shape = 21, position = position_jitterdodge(jitter.width = 1.5)) +
  scale_color_manual(values=c("grey", "grey", "grey")) +
  scale_fill_manual(values=c("white", "white", "white")) +
  labs(y=expression("Numerical DI"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
  scale_y_continuous(breaks=seq(4,12,1))


######### low and high only
ctm_pedigree_DI <- merge(family_DI, ctm_pedigree, by="AAFam")
ctm_pedigree_DI$Offspring_DI.x <- as.numeric(ctm_pedigree_DI$Offspring_DI.x)
ctm_pedigree_DI_LH <- subset(ctm_pedigree_DI, assigned_DI =="L" | assigned_DI=="H")
ctm_pedigree_DI_LH$assigned_DI <- factor(ctm_pedigree_DI_LH$assigned_DI,levels = c("L", "H"))

ctm_pedigree_DI_L <- subset(ctm_pedigree_DI, assigned_DI =="L")
mean(ctm_pedigree_DI_L$Offspring_DI.y) #6.06153
ctm_pedigree_DI_H <- subset(ctm_pedigree_DI, assigned_DI =="H")
mean(ctm_pedigree_DI_H$Offspring_DI.y) #10.96118
ctm_pedigree_DI_M <- subset(ctm_pedigree_DI, assigned_DI =="M")
mean(ctm_pedigree_DI_M$Offspring_DI.y) #8.56, range 7.1-9.95


##########
# FIG 2 plots
############
LowHigh <- subset(ctm, DI=="L" | DI=="H")
LowHigh$DI <- factor(LowHigh$DI,levels = c("L", "H"))
LowHigh$Rear_temp <- factor(LowHigh$Rear_temp,levels = c("15", "18"))

plot1 <- ggplot(data=LowHigh, aes(x=DI, y=CTM, color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  stat_summary(aes(group = Rear_temp), fun = mean, geom = "point", shape = 18, size = 4, color = "black", position = position_dodge(width = 0.75)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (B0C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))






temp_family_CTMs_DI <- merge(temp_family_CTMs, family_DI, by="AAFam")
temp_family_CTMs_DI$Offspring_DI <- as.numeric(temp_family_CTMs_DI$Offspring_DI)
temp_family_CTMs_DI_L <- subset(temp_family_CTMs_DI, Offspring_DI < 7)
temp_family_CTMs_DI_H <- subset(temp_family_CTMs_DI, Offspring_DI > 9)
temp_family_CTMs_DI_L$DI <- 'L'
temp_family_CTMs_DI_H$DI <- 'H'
temp_family_CTMs_DI_LH <- rbind(temp_family_CTMs_DI_L, temp_family_CTMs_DI_H)

temp_family_CTMs_DI_LH_Rep <- slopes_Rep %>%
  group_by(Rear_temp, AAFam) %>%
  summarise(mean_ctm = mean(mean_ctm, na.rm=TRUE))

temp_family_CTMs_DI_LH_Rep2 <- merge(temp_family_CTMs_DI_LH_Rep, ctm_pedigree, by="AAFam", all.y = F)

plot2 <- ggplot(data=temp_family_CTMs_DI_LH_Rep2, aes(x=Rear_temp.x, y=mean_ctm, group = AAFam, colour=DI)) +
  geom_point(size = 1, shape = 21) +
  geom_line(aes(color = DI)) +
  labs(y=expression("Family Mean CTM (B0C)"), x=expression("Rearing Temperature"), fill="") +
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




##################
# Models Family CTM and Family survival Correlation
##################
#load survival data from family survival script
load("smelt_survival1.RData")
load("smelt_survival2.RData")

temp_family_CTMs_Rep_big_fams$treatment_Fam <-paste(temp_family_CTMs_Rep_big_fams$AAFam, temp_family_CTMs_Rep_big_fams$Rep, temp_family_CTMs_Rep_big_fams$Rear_temp, sep = "_")
ctm.survival <- merge(temp_family_CTMs_Rep_big_fams, treatment_family_model_DI2, by="treatment_Fam")

ctm.survival.model = lm(mean_ctm ~ survival, data=ctm.survival)
Anova(ctm.survival.model, test.statistic = "F")
'''
          Sum Sq  Df F value    Pr(>F)    
survival   47.53   1  19.634 1.218e-05 ***
Residuals 951.41 393
'''
ctm.survival.model = lmer(mean_ctm ~ survival + (1|System:DI), data=ctm.survival)
Anova(ctm.survival.model, test.statistic = "F")
'''
          Sum Sq  Df F value    Pr(>F)    
              F Df Df.res    Pr(>F)    
survival 22.403  1 384.46 3.111e-06 ***
'''
plot(mean_ctm ~ survival, data=ctm.survival)
windows()
ggplot(data=ctm.survival, aes(x=survival, y=mean_ctm, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  geom_smooth(method=lm, aes(fill=Rear_temp))+
  labs(y=expression("Mean Family CTM (B0C)"), x="Percent Family Survival (%)", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
# scale_y_continuous(breaks=seq(20,30,1))


treatment_family_percent_DI$treatment_Fam <-paste(treatment_family_percent_DI$AAFam, treatment_family_percent_DI$Temp, sep = "_")
ctm.survival2 <- merge(treatment_family_percent_DI, treatment_family_model_DI2, by="treatment_Fam")

ctm.survival.model2 = lm(mean_ctm ~ survival, data=ctm.survival2)
Anova(ctm.survival.model, test.statistic = "F")
plot(mean_ctm ~ survival, data=ctm.survival)

######################## ctm slops with survival
slopes_survival <- treatment_family_percent %>%
  group_by(AAFam) %>%
  summarise(mean_survival = mean(Mean_survival, na.rm=TRUE))

slope.survival <- merge(slope_DI, slopes_survival, by="AAFam")
windows()
ggplot(data=slope.survival, aes(x=mean_survival, y=slope)) +
  geom_point(size = 1, shape = 21) +
  #scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  #scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  #geom_smooth(method=lm, aes(fill=Rear_temp))+
  labs(y=expression("Plasticity (family slope 15B0C - 18B0C)"), x="Percent Family Survival (%)", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
# scale_y_continuous(breaks=seq(20,30,1))

slope.survival <-subset(slope.survival, AAFam!="O39")
slope.survival <-subset(slope.survival, AAFam!="B4")
slope.survival.model = lm(slope ~ mean_survival, data=slope.survival)
Anova(slope.survival.model, test.statistic = "F")
'''
               Sum Sq Df F value Pr(>F)
mean_survival   0.085  1  0.0389 0.8442
Residuals     178.256 82
'''

