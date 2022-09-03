library(ggplot2)
library(reshape2)
library(jcolors)
library(ggpubr)
library(RColorBrewer)
library(plyr)
library(caret)
library(pROC)



setwd("figures/")


#### suitability all in one plot
suit <- read.csv(file="figures/suitability.csv",header=T, sep=",")
suit$variables <- factor(suit$variable, levels = c("storms","land use","SST","population","OA","DHW","overall"))
# 

### in reverse: see the number of unsuitable sites over time
suit$unsuitable <- 1- suit$suitable 
suit$pervasivee <- 1 - suit$pervasive

lag <- data.frame(x1 = c(2037, 2045), x2 = c(2048, 2073), 
                  y1 = c(0.75, 0.75), y2 = c(0.75, 0.75),
                  RCP = c("RCP 8.5 - SSP 5", "RCP 4.5 - SSP 2"),
                  variables = c("overall","overall"))

lp <- ggplot(suit, aes(x=year, y=unsuitable, group=variables))+
  geom_line(aes(color=variables,alpha=variables), size=1.5)+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7,0.7,1)) +
  geom_line(aes(y=pervasivee),linetype="dashed")


lp + facet_wrap(. ~ RCP, scales="free_x",nrow=1)+
  theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    strip.text.x = element_text(size=10, angle=0),
    strip.text.y = element_text(size=12, face="bold"),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(name="Stressor",labels=c("Storms","Land use","Population","OA","DHW","Overall"),
                     values=c("#66C2A5", "#FC8D62", "#E78AC3", "#A6D854", "#FFD92F", "#4E555A"))+
  scale_y_continuous(labels = scales::percent,expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.05,0))+
  xlab("Year")+ylab("Percent of reef sites expired")+
  scale_fill_discrete(breaks = rev(levels(suit$variables)))+
  guides(alpha=FALSE) +
  geom_segment(data = lag, aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.03, "npc"), ends = "both"),color="red")



#make barchart of percent number of expired stressors in 2100

nstress85 <- read.csv("F:/suit/rcp85-ssp5/stressor2100_count.csv",header=T)
nstress85$percent <- nstress85$nexp / sum(nstress85$nexp)
nstress85$scenario <- "RCP85-SSP5"

nstress45 <- read.csv("F:/suit/rcp45-ssp2/stressor2100_count.csv",header=T)
nstress45$percent <- nstress45$freq / sum(nstress45$freq)
nstress45$scenario <- "RCP45-SSP2"

nstress26 <- read.csv("F:/suit/rcp26-ssp1/stressor2100_count.csv",header=T)
nstress26$percent <- nstress26$freq / sum(nstress26$freq)
nstress26$scenario <- "RCP26-SSP1"


nstress <- rbind(nstress85,nstress45)
nstress <- rbind(nstress,nstress26)


lp <- ggplot(nstress, aes(x=factor(nexp2100), y=percent))+
  geom_bar(stat="identity")


lp + facet_wrap(. ~scenario,nrow=1)+
  theme_bw() + theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_blank(), 
  )+
  scale_y_continuous(labels = scales::percent_format(accuracy=1),expand = c(0,0.01)) +
  xlab("Number of stressors")+ylab("Percent of reef sites")


