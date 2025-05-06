rm(list = ls())
library(tidyverse)
library(readxl)
library(openxlsx)
library(ggpubr)

load("linelist_ref.rda")

load("single_index_caseno.rda")

load("cluster_n_case_list.rda")

load("allclustcomb_df.rda")

load("clustersize_list.rda")


# number of single index cases 
length(unique(unlist(singleindex_list)))
# 1582

# number of all cases included
allcases <- unique(c(unlist(singleindex_list),
                unlist(cluster_n_case$household$cases),
                unlist(cluster_n_case$restaurant$cases),
                unlist(cluster_n_case$office$cases),
                unlist(cluster_n_case$manual$cases),
                unlist(cluster_n_case$noso$cases),
                unlist(cluster_n_case$retail$cases),
                unlist(cluster_n_case$social$cases),
                unlist(cluster_n_case$carehome$cases)))
length(allcases)
# 8647

################################ calculate asymptomatic proportion

all_onsets <- linelist_ref %>% filter(`HK case no.` %in% allcases) %>% 
  dplyr::select(`Onset date`)
length(which(is.na(all_onsets))) # 1934
# 1934 / 8647 = 22.37%

household.clust_onsets <- linelist_ref %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$household$cases)) %>% 
  dplyr::select(`Onset date`) # 5285
length(which(is.na(household.clust_onsets))) # 1157
# 1157 / 5285         21.89%

household.single_onsets <- linelist_ref %>% 
  filter(`HK case no.` %in% singleindex_list$household) %>% 
  dplyr::select(`Onset date`) # 1536
length(which(is.na(household.single_onsets))) # 271
# 271 / 1536 = 17.64%

restaurant.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$restaurant$cases)) %>% 
  dplyr::select(`Onset date`) # 636
length(which(is.na(restaurant.clust_onsets))) # 104
# 104 / 636         16.35%

restaurant.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$restaurant) %>% 
  dplyr::select(`Onset date`) # 107
length(which(is.na(restaurant.single_onsets))) # 26
# 26 / 107 = 24.30%

office.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$office$cases)) %>% 
  dplyr::select(`Onset date`) # 235
length(which(is.na(office.clust_onsets))) # 48
# 48 / 235         20.43%

office.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$office) %>% 
  dplyr::select(`Onset date`) # 277
length(which(is.na(office.single_onsets))) # 26
# 26 / 277 = 9.39%

manual.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$manual$cases)) %>% 
  dplyr::select(`Onset date`) # 389
length(which(is.na(manual.clust_onsets))) # 105
# 105 / 389        26.99%

manual.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$manual) %>% 
  dplyr::select(`Onset date`) # 178
length(which(is.na(manual.single_onsets))) # 44
# 44 / 178 = 24.72%

noso.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$noso$cases)) %>% 
  dplyr::select(`Onset date`) # 141
length(which(is.na(noso.clust_onsets))) # 27
# 27 / 141        19.15%

noso.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$noso) %>% 
  dplyr::select(`Onset date`) # 46
length(which(is.na(noso.single_onsets))) # 2
# 2 / 46 = 4.35%


retail.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$retail$cases)) %>% 
  dplyr::select(`Onset date`) # 100
length(which(is.na(retail.clust_onsets))) # 18
# 18 / 100        18%

retail.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$retail) %>% 
  dplyr::select(`Onset date`) # 91
length(which(is.na(retail.single_onsets))) # 4
# 4 / 91 4.40%

social.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$social$cases)) %>% 
  dplyr::select(`Onset date`) # 683
length(which(is.na(social.clust_onsets))) # 196
# 196 / 683        28.70%

social.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$social) %>% 
  dplyr::select(`Onset date`) # 27
length(which(is.na(social.single_onsets))) # 4
# 4 / 27 # 14.81%

care.clust_onsets <- linelist %>% 
  filter(`HK case no.` %in% unlist(cluster_n_case$carehome$cases)) %>% 
  dplyr::select(`Onset date`) # 284
length(which(is.na(care.clust_onsets))) # 117
# 117 / 284       41.20%

care.single_onsets <- linelist %>% 
  filter(`HK case no.` %in% singleindex_list$carehome) %>% 
  dplyr::select(`Onset date`) # 36
length(which(is.na(care.single_onsets))) # 9
# 9 / 36  25%

#################################################################################

################# plot figure 2

### note this is the cluster size after checking subcluster scenarios
allclustcomb.df$category <- factor(allclustcomb.df$category,
                                   levels = c("All clusters", "Households", "Office work", "Restaurants", "Manual labour",
                                              "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes"))


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
                     axis.text = element_text(size = 7),
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


fig2b <- ggarrange(plot.fam, plot.office, plot.restaurants, plot.blue,
                   plot.retail, plot.noso, plot.indoor, plot.care,
                   nrow = 4, ncol = 2, align = "hv")
fig2b <- annotate_figure(fig2b, 
                         left = text_grob("Proportion", rot = 90),
                         bottom = text_grob("Empirical cluster size"))


pfig2 <- ggarrange(pbox, fig2b, widths = c(0.65, 0.35),
                        labels = c("a", "b"))


ggsave("Fig2_final.pdf", pfig2, width = 11.69, height = 8.27, units = "in", dpi = 600)






