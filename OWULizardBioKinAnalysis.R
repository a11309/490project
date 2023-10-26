#Load relevant libraries
library(dplyr)
library(tidyr)
library(plyr)
library(lme4)
library(car) ## ERROR
library(lmerTest)
library(ggplot2)
library(emmeans)
library(tidyverse)
library(RRPP)
library(scales)

#Clear memory
rm(list=ls(all = TRUE))

#Create theme for plots
PlotTheme <- theme_classic(base_size = 24) +
  theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) + #adjust axis and text size and position
  theme(axis.line.x=element_line(linewidth = 1.25), #Eliminates x-axis line
        axis.line.y=element_line(linewidth = 1.25), plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.ticks=element_line(size = 1.25, color="black")) + #adjust plot margins and line element size
  theme(axis.ticks.length = unit(-0.3, "cm"), legend.key.width=unit(50, "points"),
        legend.title=element_blank(),
        #legend.title=element_text(size=14, hjust=0.5),
        legend.key.height=unit(35, "points"),
        #axis.ticks.x=element_blank(), #Eliminates tick marks on x-axis
        axis.text.y=element_text(margin = margin(r = 12), color="black"),
        axis.text.x=element_text(margin = margin(t = 12), color="black"),
        #plot.title = element_text(size=16, hjust=0.5),
        plot.title = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))

####Input Data####
#Set working directory (assumes data file is available in same location as code)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
#Uploading this ridiculous file and MorphData
DLCMasterData <- read.table("DLCMasterData.csv", header=T, sep = ",", as.is=T)
MorphData <- read.table("OWULizard_Morphology.csv", header=T, sep = ",", as.is=T)
OH_Podarcis_Data <- read.table("OH_Podarcis_Data - ExpLizards_2023.csv", header=T, sep = ",", as.is=T)
#Add new columns for gravidity status sex etc. to master doc
DLCMasterDataPlus <- merge(DLCMasterData, OH_Podarcis_Data, by = "Lizard_ID", all.x = T, all.y = T)
DLCMasterDataPlus$Sex <- as.factor(DLCMasterDataPlus$Sex)
#Videos were recorded at 120 fps; Re-scale x-axis accordingly
DLCMasterDataPlus$TimeSeconds <- DLCMasterDataPlus$bodyparts_coords/120

###DLC analysis###
#Create subset for SG_PG_HL and doing the distance thing while deleting likelihood less than 0.75
# ●	PG: Pelvic Girdle
# ●	SG: Shoulder Girdle
# ●	HL: Head length; from the tip of the nose to the base of the tail

SGPGWithLikelihood <- DLCMasterDataPlus[,c(1:5,9:11,15:17,99,103:107)]
SGPGWithTwoLikelihood <- SGPGWithLikelihood[SGPGWithLikelihood$SG_likelihood > 0.75,]
SGPGWithOneLikelihood <- SGPGWithTwoLikelihood[SGPGWithTwoLikelihood$PG_likelihood > 0.75,]
SG_PG_HL <- SGPGWithOneLikelihood[SGPGWithOneLikelihood$Head_likelihood > 0.75,]
## SG_PG_HL: all 3 likelihood > 0.75

SG_PG_HL$SP_Distance <- sqrt((((SG_PG_HL$SG_x)-(SG_PG_HL$PG_x))^2)+(((SG_PG_HL$SG_y)-(SG_PG_HL$PG_y))^2))
SG_PG_HL$SH_Distance <- sqrt((((SG_PG_HL$SG_x)-(SG_PG_HL$Head_x))^2)+(((SG_PG_HL$SG_y)-(SG_PG_HL$Head_y))^2))
SG_PG_HL$PH_Distance <- sqrt((((SG_PG_HL$PG_x)-(SG_PG_HL$Head_x))^2)+(((SG_PG_HL$PG_y)-(SG_PG_HL$Head_y))^2))

#Find angle S
SG_PG_HL$Angle_S <- acos(((((SG_PG_HL$SH_Distance)^2)+((SG_PG_HL$SP_Distance)^2))-((SG_PG_HL$PH_Distance)^2))/(2*(SG_PG_HL$SH_Distance)*(SG_PG_HL$SP_Distance)))

#make only gravid females and males for now
SG_PG_HL_GFNGFM <- SG_PG_HL[,c(1,2,12:21),]
#str(SG_PG_HL_GFNGFM)
SG_PG_HL_GFNGF <- SG_PG_HL_GFNGFM[SG_PG_HL_GFNGFM$Sex != "M",] # exclude Male
#SG_PG_HL_GFNGF$Gravid
SG_PG_HL_GF <- SG_PG_HL_GFNGFM[SG_PG_HL_GFNGFM$Gravid != "N",]  # exclude nongravid
## SG_PG_HL_GF: female and gravid

#Create new column that combines sex and gravid
SG_PG_HL$Sex_Gravid <- as.factor(paste(SG_PG_HL$Sex, SG_PG_HL$Gravid, sep="_"))
str(SG_PG_HL$Sex_Gravid)

var.test(Angle_S~Sex_Gravid, data = SG_PG_HL[SG_PG_HL$Condition == "NoObs",], paired = FALSE) ## ERROR

library(car)
# Levene's test with one independent variable (three-level factor)
leveneTest(Angle_S ~ Sex_Gravid, data = SG_PG_HL[SG_PG_HL$Condition == "",]) ## ERROR

var(SG_PG_HL$Angle_S[SG_PG_HL$Sex_Gravid == "M_NA" & SG_PG_HL$Condition == "Ob"]) ## 0.0310578
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex_Gravid == "F_Y" & SG_PG_HL$Condition == "Ob"]) ## 0.04151628
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex_Gravid == "M_NA" & SG_PG_HL$Condition == "NoObs"]) ## 0.01192354
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex_Gravid == "F_Y" & SG_PG_HL$Condition == "NoObs"]) ## 0.01417495
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex_Gravid == "F_N" & SG_PG_HL$Condition == "NoObs"]) ## 0.01242974
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex_Gravid == "F_N" & SG_PG_HL$Condition == "Ob"]) ## 0.0181474

var.test(Angle_S~Sex_Gravid, data = SG_PG_HL[SG_PG_HL$Condition == "NoObs",], paired = FALSE) ## ERROR

  
#Create histogram that is fancy

# Change histogram plot fill colors by groups (NoObs)
ggplot(SG_PG_HL[SG_PG_HL$Condition == "NoObs",], aes(x=Angle_S, fill=Sex_Gravid, color=Sex_Gravid)) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_histogram(aes(y = ..density..), position="identity", alpha=0.5) +
  #geom_density(alpha=0.6, aes(y = ..scaled.., fill = Sex_Gravid)) +
  scale_fill_manual(values=c( "orange4","chartreuse4", "cornflowerblue")) +
  scale_color_manual(values=c("black", "black", "black"))+
  xlab("Angle Made by Head, Shoulder Girdle,\nand Pelvic Girdle (radians)") +
  ylab("Density") +
  #scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(limits = c(2.5,3.2)) +
  PlotTheme
  #facet_wrap(~Sex_Gravid)
ggsave(geom_histogram, file="Histogram.jpg", dpi = 250, width = 8, height = 5)

# Change histogram plot fill colors by groups (Ob)
ggplot(SG_PG_HL[SG_PG_HL$Condition == "Ob",], aes(x=Angle_S, fill=Sex_Gravid, color=Sex_Gravid)) +
  #geom_histogram(position="identity", alpha=0.5) +
  geom_histogram(aes(y = ..density..), position="identity", alpha=0.5) +
  #geom_density(alpha=0.6, aes(y = ..scaled.., fill = Sex_Gravid)) +
  scale_fill_manual(values=c( "orange4","chartreuse4", "cornflowerblue")) +
  scale_color_manual(values=c("black", "black", "black"))+
  xlab("Angle Made by Head, Shoulder Girdle,\nand Pelvic Girdle(radians)") +
  ylab("Density") +
  #scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(limits = c(2.5,3.2)) +
  PlotTheme
#facet_wrap(~Sex_Gravid)

var.test(Angle_S~Sex, data = SG_PG_HL, paired = FALSE)
var.test(Angle_S~Condition, data = SG_PG_HL, paired = FALSE)

var.test(Angle_S~Sex, data = SG_PG_HL[SG_PG_HL$Condition == "NoObs",], paired = FALSE)
str(SG_PG_HL$Sex)
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex == "M" & SG_PG_HL$Condition == "Ob"])
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex == "F" & SG_PG_HL$Condition == "Ob"])
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex == "M" & SG_PG_HL$Condition == "NoObs"])
var(SG_PG_HL$Angle_S[SG_PG_HL$Sex == "F" & SG_PG_HL$Condition == "NoObs"])

var.test(Angle_S~Sex, data = SG_PG_HL[SG_PG_HL$Condition == "Ob",], paired = FALSE)

var.test(Angle_S~Condition, data = SG_PG_HL, paired = FALSE)

#Create M_NoObs_SG_PG_HL and scatterplot for distance by frame
NoObs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Condition != "Ob",]
M_NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Sex != "F",]


#Make individual plots w/ ggplot2 for M_NoObs_SG_PG
#Create theme for plotting
plot_theme <- theme_classic(base_size = 20) +
  theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) +
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), line = element_line(linewidth = 1.25)) +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black"),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black")) +
  theme(axis.ticks.length = unit(-0.3, "cm"),
        axis.ticks=element_line(color="black")) +
  theme(panel.spacing = unit(2, units = "lines")) +
  theme(strip.background = element_blank())

M_NoObs_Angle_S_Plot <- ggplot(M_NoObs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

M_NoObs_Angle_S_Plot ## ERROR

#Create F_NoObs_SG_PG_HL and scatterplot for distance by frame
NoObs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Condition != "Ob",]
F_NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Sex != "M",]

#Make individual plots w/ ggplot2 for F_NoObs_SG_PG

F_NoObs_Angle_S_Plot <- ggplot(F_NoObs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

F_NoObs_Angle_S_Plot

#Create M_Obs_SG_PG_HL and scatterplot for distance by frame
Obs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
Obs_SG_PG_HL <- Obs_SG_PG_HL[Obs_SG_PG_HL$Condition != "NoObs",]
M_Obs_SG_PG_HL <- Obs_SG_PG_HL[Obs_SG_PG_HL$Sex != "F",]


#Make individual plots w/ ggplot2 for M_Obs_SG_PG

M_Obs_Angle_S_Plot <- ggplot(M_Obs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

M_Obs_Angle_S_Plot

#Create F_Obs_SG_PG_HL and scatterplot for distance by frame
Obs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
Obs_SG_PG_HL <- Obs_SG_PG_HL[Obs_SG_PG_HL$Condition != "NoObs",]
F_Obs_SG_PG_HL <-Obs_SG_PG_HL[Obs_SG_PG_HL$Sex != "M",]


#Make individual plots w/ ggplot2 for F_Obs_SG_PG_HL

F_Obs_Angle_S_Plot <- ggplot(F_Obs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

F_Obs_Angle_S_Plot

#Create subset for SG_PG_HL and doing the distance thing while deleting likelihood less than 0.75
SGPGTLWithLikelihood <- DLCMasterDataPlus[,c(1,2,9:11,15:17,21:23,99,103:107)]
SGPGTLWithTwoLikelihood <- SGPGTLWithLikelihood[SGPGTLWithLikelihood$SG_likelihood > 0.75,]
SGPGTLWithOneLikelihood <- SGPGTLWithTwoLikelihood[SGPGTLWithTwoLikelihood$PG_likelihood > 0.75,]
SG_PG_TL <- SGPGTLWithOneLikelihood[SGPGTLWithOneLikelihood$Tail_likelihood > 0.75,]

SG_PG_TL$SP_Distance <- sqrt((((SG_PG_TL$SG_x)-(SG_PG_TL$PG_x))^2)+(((SG_PG_TL$SG_y)-(SG_PG_TL$PG_y))^2))
SG_PG_TL$ST_Distance <- sqrt((((SG_PG_TL$SG_x)-(SG_PG_TL$Tail_x))^2)+(((SG_PG_TL$SG_y)-(SG_PG_TL$Tail_y))^2))
SG_PG_TL$PT_Distance <- sqrt((((SG_PG_TL$PG_x)-(SG_PG_TL$Tail_x))^2)+(((SG_PG_TL$PG_y)-(SG_PG_TL$Tail_y))^2))

#Find angle S
SG_PG_TL$Angle_P <- acos(((((SG_PG_TL$PT_Distance)^2)+((SG_PG_TL$SP_Distance)^2))-((SG_PG_TL$ST_Distance)^2))/(2*(SG_PG_TL$SP_Distance)*(SG_PG_TL$PT_Distance)))

hist(SG_PG_TL$Angle_P)

var.test(Angle_P~Sex, data = SG_PG_TL, paired = FALSE)
var.test(Angle_P~Condition, data = SG_PG_TL, paired = FALSE)

var.test(Angle_P~Sex, data = SG_PG_TL[SG_PG_TL$Condition == "NoObs",], paired = FALSE)
var.test(Angle_P~Sex, data = SG_PG_TL[SG_PG_TL$Condition == "Ob",], paired = FALSE)

var(SG_PG_TL$Angle_P[SG_PG_TL$Sex == "M" & SG_PG_TL$Condition == "Ob"])
var(SG_PG_TL$Angle_P[SG_PG_TL$Sex == "F" & SG_PG_TL$Condition == "Ob"])
var(SG_PG_TL$Angle_P[SG_PG_TL$Sex == "M" & SG_PG_TL$Condition == "NoObs"])
var(SG_PG_TL$Angle_P[SG_PG_TL$Sex == "F" & SG_PG_TL$Condition == "NoObs"])

var.test(Angle_S~Sex, data = SG_PG_HL[SG_PG_HL$Condition == "Ob",], paired = FALSE)

var.test(Angle_S~Condition, data = SG_PG_HL, paired = FALSE)

#Create M_NoObs_SG_PG_HL and scatterplot for distance by frame
NoObs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Condition != "Ob",]
M_NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Sex != "F",]


#Make individual plots w/ ggplot2 for M_NoObs_SG_PG
#Create theme for plotting
plot_theme <- theme_classic(base_size = 20) +
  theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) +
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), line = element_line(linewidth = 1.25)) +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black"),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black")) +
  theme(axis.ticks.length = unit(-0.3, "cm"),
        axis.ticks=element_line(color="black")) +
  theme(panel.spacing = unit(2, units = "lines")) +
  theme(strip.background = element_blank())

M_NoObs_Angle_S_Plot <- ggplot(M_NoObs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

M_NoObs_Angle_S_Plot

#Create F_NoObs_SG_PG_HL and scatterplot for distance by frame
NoObs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Condition != "Ob",]
F_NoObs_SG_PG_HL <- NoObs_SG_PG_HL[NoObs_SG_PG_HL$Sex != "M",]


#Make individual plots w/ ggplot2 for M_NoObs_SG_PG
#Create theme for plotting
plot_theme <- theme_classic(base_size = 20) +
  theme(axis.title.y=element_text(vjust=1.5), axis.title.x=element_text(vjust=0.2)) +
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), line = element_line(linewidth = 1.25)) +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) +
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black"),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black")) +
  theme(axis.ticks.length = unit(-0.3, "cm"),
        axis.ticks=element_line(color="black")) +
  theme(panel.spacing = unit(2, units = "lines")) +
  theme(strip.background = element_blank())

F_NoObs_Angle_S_Plot <- ggplot(F_NoObs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

F_NoObs_Angle_S_Plot

#Create M_Obs_SG_PG_HL and plot for distance by frame
Obs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
Obs_SG_PG_HL <- Obs_SG_PG_HL[Obs_SG_PG_HL$Condition != "NoObs",]
M_Obs_SG_PG_HL <-Obs_SG_PG_HL[Obs_SG_PG_HL$Sex != "F",]


#Make individual plots w/ ggplot2 for M_Obs_SG_PG

M_Obs_Angle_S_Plot <- ggplot(M_Obs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

M_Obs_Angle_S_Plot

#Create F_Obs_SG_PG_HL and scatterplot for distance by frame
Obs_SG_PG_HL <- SG_PG_HL[,c(1,2,12,13,17,21),]
Obs_SG_PG_HL <- Obs_SG_PG_HL[Obs_SG_PG_HL$Condition != "NoObs",]
F_Obs_SG_PG_HL <-Obs_SG_PG_HL[Obs_SG_PG_HL$Sex != "M",]


#Make individual plots w/ ggplot2 for F_Obs_SG_PG_HL

F_Obs_Angle_S_Plot <- ggplot(F_Obs_SG_PG_HL, aes(x=TimeSeconds, y=Angle_S), removePanelGrid=TRUE) + 
  geom_point(color="black", size=1.5, alpha=0.7) +
  #geom_line() +
  labs(x = "Time (s)", y = "SG-PG-HL Angle (Radians)") +
  #scale_color_manual(values = ColorPalette, name = "Color Morph") +
  #scale_fill_manual(values = FillPalette, name = "Color Morph") +
  #scale_shape_manual(values = c(21,22,23,24), name = "Color Morph") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  facet_wrap(~Lizard_ID) +
  xlim(0, 10) +
  ylim(1, 4) +
  plot_theme

F_Obs_Angle_S_Plot

#UL_SZM angle analysis
#Create subset for UL_SZM and doing the distance thing while deleting likelihood less than 0.75
UR_SZMWithLikelihood <- DLCMasterDataPlus[,c(1,2,27:35,99,103:106)]
UR_SZMWithTwoLikelihood <- UR_SZMWithLikelihood[UR_SZMWithLikelihood$UR_S_likelihood > 0.75,]
UR_SZMWithOneLikelihood <- UR_SZMWithTwoLikelihood[UR_SZMWithTwoLikelihood$UR_Z_likelihood > 0.75,]
UR_SZM <- UR_SZMWithOneLikelihood[UR_SZMWithOneLikelihood$UR_M_likelihood > 0.75,]

#Distance formula for line SZ, ZM, MS
UR_SZM$Distance_SZ <- sqrt((((UR_SZM$UR_S_x)-(UR_SZM$UR_Z_x))^2)+(((UR_SZM$UR_S_y)-(UR_SZM$UR_Z_y))^2))
UR_SZM$Distance_ZM <- sqrt((((UR_SZM$UR_M_x)-(UR_SZM$UR_Z_x))^2)+(((UR_SZM$UR_M_y)-(UR_SZM$UR_Z_y))^2))
UR_SZM$Distance_MS <- sqrt((((UR_SZM$UR_M_x)-(UR_SZM$UR_S_x))^2)+(((UR_SZM$UR_M_y)-(UR_SZM$UR_S_y))^2))
#Finding angle Z using modified Law of Cosine
UR_SZM$Angle_Z <- acos((((UR_SZM$Distance_SZ)^2)+((UR_SZM$Distance_ZM)^2)-((UR_SZM$Distance_MS)^2))/(2*(UR_SZM$Distance_SZ)*(UR_SZM$Distance_ZM)))

#Variation Test
var.test(Angle_Z~Sex, data = UR_SZM, paired = FALSE)
var.test(Angle_Z~Condition, data = UR_SZM, paired = FALSE)

var.test(Angle_Z~Sex, data = UR_SZM[UR_SZM$Condition == "NoObs",], paired = FALSE)

var(UR_SZM$Angle_Z[UR_SZM$Sex == "M" & UR_SZM$Condition == "Ob"])
var(UR_SZM$Angle_Z[UR_SZM$Sex == "F" & UR_SZM$Condition == "Ob"])
var(UR_SZM$Angle_Z[UR_SZM$Sex == "M" & UR_SZM$Condition == "NoObs"])
var(UR_SZM$Angle_Z[UR_SZM$Sex == "F" & UR_SZM$Condition == "NoObs"])

#UR_SZM angle analysis
#Create subset for UR_SZM and doing the distance thing while deleting likelihood less than 0.75
UR_SZMWithLikelihood <- DLCMasterDataPlus[,c(1,2,27:35,99,103:106)]
UR_SZMWithTwoLikelihood <- UR_SZMWithLikelihood[UR_SZMWithLikelihood$UR_S_likelihood > 0.75,]
UR_SZMWithOneLikelihood <- UR_SZMWithTwoLikelihood[UR_SZMWithTwoLikelihood$UR_Z_likelihood > 0.75,]
UR_SZM <- UR_SZMWithOneLikelihood[UR_SZMWithOneLikelihood$UR_M_likelihood > 0.75,]

#Distance formula for line SZ, ZM, MS
UR_SZM$Distance_SZ <- sqrt((((UR_SZM$UR_S_x)-(UR_SZM$UR_Z_x))^2)+(((UR_SZM$UR_S_y)-(UR_SZM$UR_Z_y))^2))
UR_SZM$Distance_ZM <- sqrt((((UR_SZM$UR_M_x)-(UR_SZM$UR_Z_x))^2)+(((UR_SZM$UR_M_y)-(UR_SZM$UR_Z_y))^2))
UR_SZM$Distance_MS <- sqrt((((UR_SZM$UR_M_x)-(UR_SZM$UR_S_x))^2)+(((UR_SZM$UR_M_y)-(UR_SZM$UR_S_y))^2))
#Finding angle Z using modified Law of Cosine
UR_SZM$Angle_Z <- acos((((UR_SZM$Distance_SZ)^2)+((UR_SZM$Distance_ZM)^2)-((UR_SZM$Distance_MS)^2))/(2*(UR_SZM$Distance_SZ)*(UR_SZM$Distance_ZM)))

#Variation Test
var.test(Angle_Z~Sex, data = UR_SZM, paired = FALSE)
var.test(Angle_Z~Condition, data = UR_SZM, paired = FALSE)

var.test(Angle_Z~Sex, data = UR_SZM[UR_SZM$Condition == "NoObs",], paired = FALSE)

var(UR_SZM$Angle_Z[UR_SZM$Sex == "M" & UR_SZM$Condition == "Ob"])
var(UR_SZM$Angle_Z[UR_SZM$Sex == "F" & UR_SZM$Condition == "Ob"])
var(UR_SZM$Angle_Z[UR_SZM$Sex == "M" & UR_SZM$Condition == "NoObs"])
var(UR_SZM$Angle_Z[UR_SZM$Sex == "F" & UR_SZM$Condition == "NoObs"])
###Clean up Morphology###
#ToDo- Replace missing toes and impute missing tail

###Morphology analysis###
t.test(TRUNK~Sex, data = MorphData, paired = FALSE)
boxplot(TRUNK~Sex, data = MorphData, xlab="Sex", ylab="Trunk Length in mm", names = c("Female", "Male"))

#Fancy boxplot for Trunk
Trunk_Boxplot <- ggplot(MorphData, aes(x=Sex, y=TRUNK, group=Sex, fill=Sex, color=Sex), horizontal=FALSE, removePanelGrid=TRUE) + 
  geom_boxplot(linewidth = 1, alpha = 0.75, color="black", width = 0.5, outlier.shape = NA) +
  #geom_line(aes(group=Lizard_ID, alpha=0.5), color="black",  linewidth = 1, position = position_dodge(0.2)) + 
  #geom_point(aes(alpha=1), size = 5, position = position_dodge(0.2)) +
  #geom_jitter(position=position_jitter(0.1)) +
  scale_fill_manual(values=c("chartreuse4", "cornflowerblue")) +
  scale_color_manual(values=c("chartreuse4", "cornflowerblue")) +
  geom_boxplot(size = 1, alpha = 0, color="black", width = 0.5, outlier.shape = NA) +
  #stat_boxplot(size = 2, color="black", geom = "errorbar", width = 0.3) +
  labs(x = element_blank()) +
  scale_x_discrete(labels=c("Female", "Male")) +
  ylab("Trunk Length (mm)") +
  ylim(20, 40) +
  #geom_hline(yintercept = 0.1, linewidth=2, linetype="dashed", color = "orange") +
  #geom_hline(yintercept = 0.8, linewidth=2, linetype="dashed", color = "red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color="black", size= 24),
        axis.text.y = element_text(color="black", size = 22),
        axis.title.y=element_text(vjust=1.5, size = 24),
        axis.ticks=element_line(color="black"), axis.ticks.length = unit(-0.3, "cm"),
        axis.title=element_text(size=18), legend.position="none")
Trunk_Boxplot


