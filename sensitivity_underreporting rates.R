rm(list = ls())
library(ggplot2)
library(dplyr)
library(scales)
source("tupper.R")
load("size_by_setting_0410.rda")
clustall <- read.csv("allclustersize.csv")

res.all <- optim_fit_model(init = c(0.5, 0.5),
                           inputdata = list(q1 = 0.8,
                                            q2 = 0.7,
                                            clusters = clustall$clustsize))

outputres <- function(modelres){
  data.frame(
    dispersion_k = modelres$a_n_b.est[1],
    meanclust = 1 + modelres$a_n_b.est[1] * modelres$a_n_b.est[2],
    Rc = modelres$a_n_b.est[1] * modelres$a_n_b.est[2]
  )
}

outputres(res.all)

est_k_allsetting <- outputres(res.all)[1]
est_R_allsetting <- outputres(res.all)[3]

q1_opt <- seq(0.5, 1.0, 0.1)
q2_opt <- seq(0.1, 1.0, 0.1)

est_k_opt <- est_R_opt <- matrix(nrow = 6, ncol = 10, byrow = T)

for(i in 1:6){
  for(j in 1:10){
    fittmp <- optim_fit_model(init = c(0.5, 0.5),
                                         inputdata = list(q1 = q1_opt[i],
                                                          q2 = q2_opt[j],
                                                          clusters = clustall$clustsize))
    outtmp <- outputres(fittmp)
    est_k_opt[i, j] <- as.numeric(outtmp[1])
    est_R_opt[i, j] <- as.numeric(outtmp[3])
    }
}

data <- expand.grid(X=q1_opt, Y=q2_opt)
data$Z <- as.numeric(est_k_opt)
data$Z2 <- as.numeric(est_R_opt)

# Give extreme colors:
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() + xlab("q1 (probability of being reported)") + ylab("q2 (probability of being related to clusters)") +
  scale_fill_gradient(low="white", high="blue", name = "Dispersion parameter k") + theme(legend.position ="top")

ggplot(data, aes(X, Y, fill= Z2)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") 


# build a function
possible_para <- function(clustersize){
  q1_opt <- seq(0.5, 1.0, 0.1)
  q2_opt <- seq(0.1, 1.0, 0.1)
  est_k_opt <- est_R_opt <- matrix(nrow = 6, ncol = 10, byrow = T)
  for(i in 1:6){
    for(j in 1:10){
      fittmp <- optim_fit_model(init = c(0.5, 0.5),
                                inputdata = list(q1 = q1_opt[i],
                                                 q2 = q2_opt[j],
                                                 clusters = clustersize))
      outtmp <- outputres(fittmp)
      est_k_opt[i, j] <- as.numeric(outtmp[1])
      est_R_opt[i, j] <- as.numeric(outtmp[3])
    }
  }
  return(list(
    k_values = est_k_opt,
    R_values = est_R_opt
  ))
}

clust_size_list <- list(
  allclusters = clustall$clustsize,
  household = size_by_setting$household$clustsize,
  carehome = size_by_setting$carehome$clustsize,
  restaurant = size_by_setting$restaurant$clustsize,
  noso = size_by_setting$nosocomial$clustsize,
  blue = size_by_setting$blue$clustsize,
  office = size_by_setting$office$clustsize,
  GCLA = size_by_setting$retail$clustsize,
  ICFS = size_by_setting$indoor$clustsize
)


progressbar <- txtProgressBar(min = 0, max = 9, style = 3)
tstart <- Sys.time()
k_R_list <- vector("list", 9)
for(m in 1:9){
  k_R_list[[m]] <- possible_para(clust_size_list[[m]])
  setTxtProgressBar(progressbar, m)
}
tend <- Sys.time()
tend - tstart

save(k_R_list, file = "sensitivity_k_R_20230815.rda")
save(k_R_list, file = "sensitivity_k_R_20240411.rda")

# load("sensitivity_k_R_20230815.rda")
load("sensitivity_k_R_20240411.rda")

plot_k_list <- vector("list", 9)
plot_R_list <- vector("list", 9)
plot_titles <- c("All clusters", "Households", "Care homes", "Restaurants",
                 "Nosocomial", "Blue-collar", "Office work", "Retail & leisure", "Close-social indoor")

for(i in 1:9){
  data <- expand.grid(X=q1_opt, Y=q2_opt)
  data$Z <- as.numeric(k_R_list[[i]][[1]])
  data$Z2 <- as.numeric(k_R_list[[i]][[2]])
  plot_k_list[[i]] <- ggplot(data, aes(X, Y, fill= Z)) + 
    geom_tile() + xlab("q1 (probability of being reported)") + ylab("q2 (probability of being related to clusters)") +
    scale_fill_gradient(low="white", high="blue", name = "Dispersion parameter k", n.breaks = 4) + 
    theme(legend.position ="top", legend.text=element_text(size=6),legend.title = element_text(size=8)) +
    scale_x_continuous(breaks = seq(0.5, 1.0, 0.1)) + scale_y_continuous(breaks = seq(0.1, 1.0, 0.1)) + 
    theme_bw() + theme(panel.border = element_blank(),
                       axis.line = element_line(),
                       legend.position = "top") +
    ggtitle(plot_titles[i])
  plot_R_list[[i]] <- ggplot(data, aes(X, Y, fill= Z2)) + 
    geom_tile() + xlab("q1 (probability of being reported)") + ylab("q2 (probability of being related to clusters)") +
    scale_fill_gradient(low="white", high="red", name = "Expected no. new infections", n.breaks = 4) + 
    theme(legend.position ="top", legend.text=element_text(size=6),legend.title = element_text(size=8)) +
    scale_x_continuous(breaks = seq(0.5, 1.0, 0.1)) + scale_y_continuous(breaks = seq(0.1, 1.0, 0.1)) + 
    theme_bw() + theme(panel.border = element_blank(),
                       axis.line = element_line(),
                       legend.position = "top") +
    ggtitle(plot_titles[i])
}

save(plot_k_list, file = "sensitivity_k_plot_20230815.rda")
save(plot_R_list, file = "sensitivity_R_plot_20230815.rda")


save(plot_k_list, file = "sensitivity_k_plot_20240411.rda")
save(plot_R_list, file = "sensitivity_R_plot_20240411.rda")


library(ggpubr)

p1 <- ggarrange(plot_k_list[[1]], plot_R_list[[1]], nrow = 1, ncol = 2)
p2 <- ggarrange(plot_k_list[[2]], plot_R_list[[2]], nrow = 1, ncol = 2)
p3 <- ggarrange(plot_k_list[[3]], plot_R_list[[3]], nrow = 1, ncol = 2)
p4 <- ggarrange(plot_k_list[[4]], plot_R_list[[4]], nrow = 1, ncol = 2)
p5 <- ggarrange(plot_k_list[[5]], plot_R_list[[5]], nrow = 1, ncol = 2)
p6 <- ggarrange(plot_k_list[[6]], plot_R_list[[6]], nrow = 1, ncol = 2)
p7 <- ggarrange(plot_k_list[[7]], plot_R_list[[7]], nrow = 1, ncol = 2)
p8 <- ggarrange(plot_k_list[[8]], plot_R_list[[8]], nrow = 1, ncol = 2)
p9 <- ggarrange(plot_k_list[[9]], plot_R_list[[9]], nrow = 1, ncol = 2)

p_heatmap_k <- ggarrange(plot_k_list[[1]], plot_k_list[[2]], plot_k_list[[7]],
                         plot_k_list[[4]], plot_k_list[[6]], plot_k_list[[8]],
                         plot_k_list[[5]], plot_k_list[[9]], plot_k_list[[3]],
                       nrow = 3, ncol = 3, labels = letters[1:9])
p_heatmap_k
ggsave("./event cluster analysis/manuscript/FigS1_20240116.pdf",p_heatmap_k, width = 16, height = 12, units = "in")
ggsave("./event cluster analysis/manuscript/FigS1_20240411.pdf",p_heatmap_k, width = 16, height = 12, units = "in")
ggsave("./event cluster analysis/manuscript/FigS1_20240411.pdf",p_heatmap_k, width = 16, height = 12, units = "in")


p_heatmap_R <- ggarrange(plot_R_list[[1]], plot_R_list[[2]], plot_R_list[[7]],
                         plot_R_list[[4]], plot_R_list[[6]], plot_R_list[[8]],
                         plot_R_list[[5]], plot_R_list[[9]], plot_R_list[[3]],
                         nrow = 3, ncol = 3, labels = letters[1:9])
p_heatmap_R
ggsave("./event cluster analysis/manuscript/FigS2_20240116.pdf",p_heatmap_R, width = 16, height = 12, units = "in")
ggsave("./event cluster analysis/manuscript/FigS2_20240321.pdf",p_heatmap_R, width = 16, height = 12, units = "in")
ggsave("./event cluster analysis/manuscript/FigS2_20240411.pdf",p_heatmap_R, width = 16, height = 12, units = "in")


