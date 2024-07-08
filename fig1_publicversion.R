rm(list = ls())

library(ggplot2)
library(dplyr)
library(scales)
load("size_by_setting_0410.rda")
clustall <- read.csv("allclustersize.csv")

## allclust size combine
allclustcomb.df <- rbind(clustall[,-1], 
                         size_by_setting$household,
                         size_by_setting$carehome,
                         size_by_setting$restaurant,
                         size_by_setting$nosocomial,
                         size_by_setting$blue,
                         size_by_setting$office,
                         size_by_setting$indoor,
                         size_by_setting$retail)


allclustcomb.df$category <- factor(allclustcomb.df$category,
                                   levels = c("All clusters", "Households", "Office work", "Restaurants", "Manual labour",
                                              "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes"))


pbox <- ggplot(data = allclustcomb.df, aes(y = clustsize,  x = category)) + 
  geom_boxplot(aes(color = category, fill = category), alpha = 0.25) + 
  scale_y_log10(breaks = c(1, 2, 3, 5, 10, 50, 100, 300, 500)) +
  geom_point(aes(y = clustsize,  x = category, color = category), position = "jitter", alpha = 0.1) +
  scale_color_manual(values=c("#000000", "#32CD32", 
                              "#808080", "#FF8C00",
                              "#1E90FF", "#008080",
                              "#800080", "#FF1493", 
                              "#8B4513"))+
  scale_fill_manual(values=c("#000000", "#32CD32", 
                             "#808080", "#FF8C00",
                             "#1E90FF", "#008080",
                             "#800080", "#FF1493", 
                             "#8B4513"))+
  scale_x_discrete(labels = c("All clusters\n2214 + 1576",
                              "Households\n1776 + 1540", 
                              "Office work\n76 + 270", 
                              "Restaurants\n177 + 105", 
                              "Manual labour\n77 + 157",
                              "Retail & leisure\n18 + 66", 
                              "Nosocomial\n34 + 45", 
                              "Close-social indoor\n33 + 31", "Care homes\n23 + 34"))+
  xlab("Cluster setting (count)\nNo. of clusters involved at least 2 cases + No. of single index cases") + 
  ylab("Empirical cluster size") + 
  theme_bw() + theme(legend.position = "none",
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10)) 

pbox

### frequency by setting
unique(size_by_setting$household$clustsize)

## create a function to provide plots
plotfun <- function(data, coloropt, clustname){
  dftmp <- data %>% 
    mutate(sizecate = case_when(
      clustsize == 1 ~ "1",clustsize == 2 ~ "2",clustsize == 3 ~ "3",
      clustsize == 4 ~ "4",clustsize == 5 ~ "5",clustsize == 6 ~ "6",
      clustsize == 7 ~ "7",clustsize == 8 ~ "8",clustsize == 9 ~ "9",
      clustsize >= 10 ~ "10+",
    )
    )
  dt <- table(dftmp$sizecate) %>% as.data.frame() %>% 
    mutate(prop = Freq/sum(Freq)) %>% mutate(
      xmin = case_when(
        Var1 == "1" ~ 1 - 0.35,Var1 == "2" ~ 2 - 0.35,Var1 == "3" ~ 3 - 0.35,
        Var1 == "4" ~ 4 - 0.35,Var1 == "5" ~ 5 - 0.35,Var1 == "6" ~ 6 - 0.35,
        Var1 == "7" ~ 7 - 0.35,Var1 == "8" ~ 8 - 0.35,Var1 == "9" ~ 9 - 0.35,
        Var1 == "10+" ~ 10 - 0.35,
      ),
      xmax = case_when(
        Var1 == "1" ~ 1 + 0.35,Var1 == "2" ~ 2 + 0.35,Var1 == "3" ~ 3 + 0.35,
        Var1 == "4" ~ 4 + 0.35,Var1 == "5" ~ 5 + 0.35,Var1 == "6" ~ 6 + 0.35,
        Var1 == "7" ~ 7 + 0.35,Var1 == "8" ~ 8 + 0.35,Var1 == "9" ~ 9 + 0.35,
        Var1 == "10+" ~ 10 + 0.35,
      )
    )
  plot <- ggplot(dt) +
    geom_rect(aes(xmin = xmin,  xmax = xmax, ymin = 0, ymax = prop),
              fill = coloropt, alpha = 0.5) +
    # xlab("Empirical cluster size") + ylab("Proportion") +
    scale_x_continuous(breaks = seq(1, 10), labels = c(seq(1, 9), "10+")) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) + 
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid = element_blank(),
                       axis.line = element_line()) +
    ggtitle(clustname)
  
  return(plot)
}


plot.fam <- plotfun(size_by_setting$household,
                    coloropt = "#32CD32",
                    clustname = "Households")

plot.office <-  plotfun(size_by_setting$office,
                        coloropt = "#808080",
                        clustname = "Office work")

plot.restaurants <-  plotfun(size_by_setting$restaurant,
                             coloropt = "#FF8C00",
                             clustname = "Restaurants")

plot.blue <- plotfun(size_by_setting$blue,
                     coloropt = "#1E90FF",
                     clustname = "Manual labour")

plot.retail <- plotfun(size_by_setting$retail,
                       coloropt = "#008080",
                       clustname = "Retail & leisure")

plot.noso <- plotfun(size_by_setting$nosocomial,
                     coloropt = "#800080",
                     clustname = "Nosocomial")

plot.indoor <- plotfun(size_by_setting$indoor,
                       coloropt = "#FF1493",
                       clustname = "Close-social indoor")

plot.care <- plotfun(size_by_setting$carehome,
                     coloropt = "#8B4513",
                     clustname = "Care homes")

fig1b <- ggarrange(plot.fam, plot.office, plot.restaurants, plot.blue,
                   plot.retail, plot.noso, plot.indoor, plot.care,
                   nrow = 4, ncol = 2, align = "hv")
fig1b <- annotate_figure(fig1b, 
                         left = text_grob("Proportion", rot = 90),
                         bottom = text_grob("Empirical cluster size"))


pfig1_main <- ggarrange(pbox, fig1b, widths = c(0.6, 0.4),
                        labels = c("a", "b"))


# ggsave("Fig1_main.pdf", pfig1_main, width = 13.5, height = 7.5, units = "in", dpi = 900)

