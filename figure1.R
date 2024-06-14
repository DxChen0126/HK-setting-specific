
library(ggplot2)
library(dplyr)
library(scales)
source("tupper.R")
load("allclustsize.rda")

## household size
fam.clust <- read.csv("fam_event_freq.csv")
fam.single <- read.csv("fam_single_event_freq.csv")

famsize <- c(fam.clust$clust_freq, rep(1, nrow(fam.single)))
head(famsize)

famsize.df <- data.frame(
 clustsize = famsize,
 category = rep("Households", length(famsize))
)

fam.single.caseno <- fam.single$HK.case.no.


## care home
care.clust <- read.csv("care_event_freq.csv")
care.single <- read.csv("care_single_event_freq.csv")
which(care.single$assumed.event == "white")
caresize <- c(care.clust$clust_freq, rep(1, nrow(care.single) - 1))
head(caresize)

caresize.df <- data.frame(
  clustsize = caresize,
  category = rep("Care homes", length(caresize))
)

care.single.caseno <- care.single$HK.case.no.[-which(care.single$assumed.event == "white")]


## nosocomial; blue collar; office work; catering; they need reclassification
## and had to read all of them then distinguish each other
noso.clust <- read.csv("noso_event_freq.csv")
noso.single <- read.csv("noso_single_event_freq.csv")

blue.clust <- read.csv("blue_event_freq.csv")
blue.single <- read.csv("blue_single_event_freq.csv")

white.clust <- read.csv("white_event_freq.csv")
white.single <- read.csv("white_single_event_freq.csv")

cater.clust <- read.csv("cater_event_freq.csv")
cater.single <- read.csv("cater_single_event_freq.csv")

dine.clust <- read.csv("dine_event_freq.csv")
dine.single <- read.csv("dine_single_event_freq.csv")

unique(cater.single$reclassify)

dinesize <- c(dine.clust$clust_freq, rep(1, nrow(dine.single)),
              rep(1, length(which(cater.single$reclassify=="DINE"))))
head(dinesize)

dine.single.caseno <- c(dine.single$HK.case.no., 
                        cater.single$HK.case.no.[which(cater.single$reclassify=="DINE")])


## restaurant
dinesize.df <- data.frame(
  clustsize = dinesize,
  category = rep("Restaurants", length(dinesize))
)

## Noso
nososize <- c(noso.clust$clust_freq, rep(1, nrow(noso.single)),
              rep(1, length(which(blue.single$reclassify == "noso"))))
head(nososize)

nososize.df <- data.frame(
  clustsize = nososize,
  category = rep("Nosocomial", length(nososize))
)

noso.single.caseno <- c(noso.single$HK.case.no.,
                        blue.single$HK.case.no.[which(blue.single$reclassify == "noso")])


## blue
unique(blue.clust$reclassify)
ind.blue1 <- which(blue.clust$reclassify=="")
unique(blue.single$reclassify)
ind.blue2 <- which(blue.single$reclassify=="")
unique(white.clust$reclassify)
ind.blue3 <- which(white.clust$reclassify=="BW")
unique(white.single$reclassify)
ind.blue4 <- which(white.single$reclassify=="BW")
bluesize <- c(blue.clust[ind.blue1,]$clust_freq,
              white.clust[ind.blue3,]$clust_freq,
              rep(1, length(ind.blue2)), rep(1, length(ind.blue4)))

bluesize.df <- data.frame(
  clustsize = bluesize,
  category = rep("Blue-collar", length(bluesize))
)

blue.single.caseno <- c(blue.single$HK.case.no.[ind.blue2],
                        white.single$HK.case.no.[ind.blue4])


## Office work
unique(white.clust$reclassify)
ind.white1 <- which(white.clust$reclassify=="")
unique(white.single$reclassify)
ind.white2 <- which(white.single$reclassify=="")
unique(cater.single$reclassify)
ind.white3 <- which(cater.single$reclassify=="OW")

whitesize <- c(white.clust[ind.white1,]$clust_freq,
               rep(1, length(ind.white2)), rep(1, length(ind.white3)))
whitesize.df <- data.frame(
  clustsize = whitesize,
  category = rep("Office work", length(whitesize))
)

white.single.caseno <- c(white.single$HK.case.no.[ind.white2],
                         cater.single$HK.case.no.[ind.white3])


## indoor catering
ind.ICFS1 <- which(cater.clust$reclassify == "ICFS")
ind.ICFS2 <- which(cater.single$reclassify == "ICFS")
ICFSsize <- c(cater.clust[ind.ICFS1,]$clust_freq,
              rep(1, length(ind.ICFS2)))
ICFSsize.df <- data.frame(
  clustsize = ICFSsize,
  category = rep("Close-social indoor", length(ICFSsize))
)

ICFS.single.caseno <- cater.single$HK.case.no.[ind.ICFS2]

## retail (general consumption)
unique(cater.clust$reclassify)
unique(cater.single$reclassify)
unique(blue.clust$reclassify)
ind.GCLA1 <- which(cater.clust$reclassify == "GCLA")
ind.GCLA2 <- which(cater.single$reclassify == "GCLA")
ind.GCLA3 <- which(blue.clust$reclassify == "GCLA")
GCLAsize <- c(cater.clust[ind.GCLA1,]$clust_freq,
              blue.clust[ind.GCLA3,]$clust_freq,
              rep(1, length(ind.GCLA2)))
GCLAsize.df <- data.frame(
  clustsize = GCLAsize,
  category = rep("Retail & leisure", length(GCLAsize))
)

GCLA.single.caseno <- cater.single$HK.case.no.[ind.GCLA2]
                        

allsingle.caseno <- c(fam.single.caseno,
                      dine.single.caseno,
                      noso.single.caseno,
                      blue.single.caseno,
                      white.single.caseno,
                      care.single.caseno,
                      ICFS.single.caseno,
                      GCLA.single.caseno)

length(unique(allsingle.caseno)) # 1577
library(tidyverse)

tb <- as.data.frame(table(allsingle.caseno))
head(tb)
ind.dup <- which(tb$Freq > 1)
unique(tb[ind.dup,]$Freq)

1623 - 1577

allsize.df <- data.frame(
  clustsize = clustsize.incsingle,
  category = rep("All clusters", length(clustsize.incsingle))
)

## allclust size combine
allclustcomb.df <- rbind(allsize.df, famsize.df, caresize.df, dinesize.df, nososize.df,
                         bluesize.df, whitesize.df, ICFSsize.df, GCLAsize.df)


nrow(allsize.df)
nrow(famsize.df)


allclustcomb.df$category <- factor(allclustcomb.df$category,
                                   levels = c("All clusters", "Households", "Office work", "Restaurants", "Blue-collar",
                                              "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes"))


pbox <- ggplot(data = allclustcomb.df, aes(y = clustsize,  x = category)) + 
  geom_boxplot(aes(color = category)) + 
  scale_y_log10(breaks = c(1, 5, 10, 50, 100, 300)) +
  geom_point(aes(y = clustsize,  x = category, color = category), position = "jitter", alpha = 0.1) +
  scale_color_manual(values=c("#000000", "#32CD32", 
                              "#808080", "#FF8C00",
                              "#1E90FF", "#008080",
                              "#800080", "#FF1493", 
                              "#8B4513"))+
  scale_x_discrete(labels = c("All clusters\n4264",
                              "Households\n 3327", 
                              "Office work\n346", 
                              "Restaurants\n283", 
                              "Blue-collar\n234",
                              "Retail & leisure\n84", 
                              "Nosocomial\n79", 
                              "Close-social indoor\n64", "Care homes\n46"))+
  xlab(" ") + ylab("Cluster size") + 
  theme_bw() + theme(legend.position = "none",
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10)) 


clustcountdf <- data.frame(
  freq = c(nrow(allsize.df), nrow(famsize.df),  nrow(whitesize.df), nrow(dinesize.df), nrow(bluesize.df), 
           nrow(GCLAsize.df), nrow(nososize.df), nrow(ICFSsize.df), nrow(caresize.df)),
  category = c("All clusters","Households", "Office work", "Restaurants", "Blue-collar",
               "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes"),
  Xind = seq(1, 9, 1)
) 

# clustcountdf$category <- factor(clustcountdf$category,
#                                    levels = c("All clusters","Households", "Office work", "Restaurants", "Blue-collar",
#                                               "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes"))
# 
# 
# pcountbar <- ggplot(data = clustcountdf) +
#   geom_rect(aes(xmin = Xind - 0.25, xmax = Xind + 0.25, ymin = 0, ymax = freq, fill = category), alpha = 0.75) +
#   scale_y_log10() + scale_fill_manual(values=c("#000000", "#32CD32", 
#                                                "#808080", "#FF8C00",
#                                                "#1E90FF", "#008080",
#                                                "#800080", "#FF1493", 
#                                                "#8B4513"))+
#   xlab("Cluster category") + ylab("Cluster frequency") + 
#   scale_x_continuous(breaks = seq(1, 9), labels = c("All clusters","Households", "Office work", "Restaurants", "Blue-collar",
#                                                     "Retail &\nleisure", "Nosocomial", "Close-social\nindoor", "Care homes"))+
#   theme_bw() + theme(legend.position = "none",
#                      panel.grid.minor = element_blank(),
#                      axis.text = element_text(size = 8),
#                      axis.title = element_text(size = 10)) 
# 
# 
# 
# 
# library(ggpubr)

pfig1_main <- ggarrange(pviolin, pcountbar, nrow = 2, ncol = 1, labels = letters[1:2], align = "v")

ggsave("./event cluster analysis/manuscript/Fig1_20240116.pdf", pfig1_main, width = 8, height = 7.5, units = "in", dpi = 900)

ggsave("./event cluster analysis/manuscript/Fig1_20240122.pdf", pviolin, width = 9, height = 6, units = "in", dpi = 900)

ggsave("./event cluster analysis/manuscript/Fig1_20240321.pdf", pbox, width = 9, height = 6, units = "in", dpi = 900)





















