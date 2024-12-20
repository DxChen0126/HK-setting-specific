library(tidyverse)
library(readxl)
library(openxlsx)
library(ggpubr)


### these two are the anonymized data contained cluster name, cases involved and 
### cluster size
load("cluster_n_case_list.rda")
load("clustersize_list.rda")

### this is the processed dataframe for visualization and model estimation
load("allclustcomb_df.rda")

pbox <- ggplot(data = allclustcomb.df, aes(y = size,  x = category)) + 
  geom_boxplot(aes(color = category, fill = category), alpha = 0.25) + 
  scale_y_log10(breaks = c(1, 2, 3, 5, 10, 50, 100, 300, 500)) +
  geom_point(aes(y = size,  x = category, color = category),alpha = 0.1, position = "jitter") +
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
  scale_x_discrete(labels = c("All clusters\n2228 + 1582",
                              "Households\n1782 + 1536", 
                              "Office work\n88 + 277", 
                              "Restaurants\n175 + 107", 
                              "Manual labour\n75 + 178",
                              "Retail & leisure\n17 + 91", 
                              "Nosocomial\n34 + 46", 
                              "Close-social indoor\n34 + 27", "Care homes\n23 + 36"))+
  xlab("No. of clusters involved at least 2 cases + No. of single index cases
       under each specified transmission setting") + 
  ylab("Empirical cluster size") + 
  theme_bw() + theme(legend.position = "none",
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size = 8),
                     axis.title = element_text(size = 10)) 

pbox

## create a function to provide plots
plotfun <- function(data, coloropt, clustname){
  dftmp <- data %>% 
    mutate(sizecate = case_when(
      clustersize == 1 ~ "1",clustersize == 2 ~ "2",clustersize == 3 ~ "3",
      clustersize == 4 ~ "4",clustersize == 5 ~ "5",clustersize == 6 ~ "6",
      clustersize == 7 ~ "7",clustersize == 8 ~ "8",clustersize == 9 ~ "9",
      clustersize >= 10 ~ "10+",
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
    scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2)) + 
    theme_bw() + theme(panel.border = element_blank(),
                       panel.grid = element_blank(),
                       axis.line = element_line()) +
    ggtitle(clustname)
  
  return(plot)
}


plot.fam <- plotfun(clustersize_list$household,
                    coloropt = "#32CD32",
                    clustname = "Households")

plot.office <-  plotfun(clustersize_list$office,
                        coloropt = "#808080",
                        clustname = "Office work")

plot.restaurants <-  plotfun(clustersize_list$restaurant,
                             coloropt = "#FF8C00",
                             clustname = "Restaurants")

plot.blue <- plotfun(clustersize_list$manual_labour,
                     coloropt = "#1E90FF",
                     clustname = "Manual labour")

plot.retail <- plotfun(clustersize_list$retail,
                       coloropt = "#008080",
                       clustname = "Retail & leisure")

plot.noso <- plotfun(clustersize_list$nosocomial,
                     coloropt = "#800080",
                     clustname = "Nosocomial")

plot.indoor <- plotfun(clustersize_list$social,
                       coloropt = "#FF1493",
                       clustname = "Close-social indoor")

plot.care <- plotfun(clustersize_list$carehome,
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


ggsave("Fig1_main_20241212.pdf", pfig1_main, width = 13.5, height = 7.5, units = "in", dpi = 900)
