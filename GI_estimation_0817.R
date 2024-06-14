rm(list = ls())

library(dplyr)
library(foreach)
library(doParallel)

ncores <- detectCores()

load("./sample infection record/main_GI_100times_sample_20230803.rda")

source("mixture model 0801.R")


progressbar <- txtProgressBar(min = 0, max = 100, style = 3)
tstart <- Sys.time()
SampleGT <- sampleGT_record[[8]]
GTest <- vector("list", 100)
for(j in 1:100){
  set.seed(j)
  tempdata <- SampleGT[[1 + (j-1)*2]]
  tempest <- GI_mix_est_revise(data = tempdata, N = 100, startmu = 10, startsig = 10)
  save(tempest, file = paste0("./sample infection record/GCLA/GI_est_GCLA_", j, "th_est.rda"))
  # do parallel for the bootstap part
  myCluster <- makeCluster(ncores - 1, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  estboot <- foreach(i = 1:100, .combine = 'append',
                     .packages=c("dplyr"), .errorhandling = 'pass') %dopar% {
                       GI_mix_est_revise(data = sample(tempdata, length(tempdata), replace = T),
                                         N = 100, startmu = 10, startsig = 10)
                     }
  save(estboot, file = paste0("./sample infection record/GCLA/GI_boot_GCLA_", j, "th_sample.rda"))
  # failed at the 31th bootstrap, will amend later
  stopCluster(myCluster)
  print(j)
  GTest[[j]] <- list(estimate = tempest,
                     bootstap.res = estboot)
  setTxtProgressBar(progressbar, j)
}

tend <- Sys.time()
tend - tstart

save(GTest, file = "./sample infection record/main_GI_est_GCLA_20230817.rda")


### apply rubin's rule on log of the estimates?

### still makes no sense

### a safer way might be get the 95% credible interval from all the bootstrap record

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)
tstart <- Sys.time()
SampleGT <- sampleGT_record[[7]]
GTest <- vector("list", 100)
for(j in 1:100){
  set.seed(j)
  tempdata <- SampleGT[[1 + (j-1)*2]]
  tempest <- GI_mix_est_revise(data = tempdata, N = 100, startmu = 10, startsig = 10)
  save(tempest, file = paste0("./sample infection record/ICFS/GI_est_ICFS_", j, "th_est.rda"))
  # do parallel for the bootstrap part
  myCluster <- makeCluster(ncores - 1, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  estboot <- foreach(i = 1:100, .combine = 'append',
                     .packages=c("dplyr"), .errorhandling = 'pass') %dopar% {
                       GI_mix_est_revise(data = sample(tempdata, length(tempdata), replace = T),
                                         N = 100, startmu = 10, startsig = 10)
                     }
  save(estboot, file = paste0("./sample infection record/ICFS/GI_boot_ICFS_", j, "th_sample.rda"))
  # failed at the 31th bootstrap, will amend later
  stopCluster(myCluster)
  print(j)
  GTest[[j]] <- list(estimate = tempest,
                     bootstap.res = estboot)
  setTxtProgressBar(progressbar, j)
}

tend <- Sys.time()
tend - tstart

save(GTest, file = "./sample infection record/main_GI_est_ICFS_20230822.rda")


####### bootstrap CI

## 1) ICFS
obtain_est <- function(estdata){
  mu.est <- sigma.est <- w1.est <- w2.est <- w3.est <- numeric(100)
  boot_mu.est <- boot_sigma.est <- boot_w1.est <- boot_w2.est <- boot_w3.est <- vector("list", 100)
  for(i in 1:100){
    tmp <- estdata[[i]]
    mu.est[i] <- tmp[[1]][1]
    sigma.est[i] <- tmp[[1]][2]
    w1.est[i] <- tmp[[1]][3]
    w2.est[i] <- tmp[[1]][4]
    w3.est[i] <- tmp[[1]][5]
    boot_mu.est[[i]] <- tmp[[2]][seq(1, 500, 5)]
    boot_sigma.est[[i]] <- tmp[[2]][seq(2, 500, 5)]
    boot_w1.est[[i]] <- tmp[[2]][seq(3, 500, 5)]
    boot_w2.est[[i]] <- tmp[[2]][seq(4, 500, 5)]
    boot_w3.est[[i]] <- tmp[[2]][seq(5, 500, 5)]
  }
  return(data.frame(mu.pool = mean(mu.est),
                    mu.lb = quantile(unlist(boot_mu.est), 0.025),
                    mu.ub = quantile(unlist(boot_mu.est), 0.975),
                    sigma.pool = mean(sigma.est),
                    sigma.lb = quantile(unlist(boot_sigma.est), 0.025),
                    sigma.ub = quantile(unlist(boot_sigma.est), 0.975),
                    w1.pool = mean(w1.est),
                    w1.lb = quantile(unlist(boot_w1.est), 0.025),
                    w1.ub = quantile(unlist(boot_w1.est), 0.975),
                    w2.pool = mean(w2.est),
                    w2.lb = quantile(unlist(boot_w2.est), 0.025),
                    w2.ub = quantile(unlist(boot_w2.est), 0.975),
                    w3.pool = mean(w3.est),
                    w3.lb = quantile(unlist(boot_w3.est), 0.025),
                    w3.ub = quantile(unlist(boot_w3.est), 0.975)
  ))
  
}

load("./sample infection record/main_GI_est_ICFS_20230822.rda")
pool_ICFS <- obtain_est(GTest)

# 2) GCLA
load("./sample infection record/main_GI_est_GCLA_20230817.rda")
pool_GCLA <- obtain_est(GTest)

# 3) care homes
load("Main_GI_est_carehome_20230817.rda")
pool_carehome <- obtain_est(GTest)

# 4) restaurants
load("Main_GI_est_dining_20230818.rda")
pool_restaurant <- obtain_est(GTest)

# 5) noso
load("Main_GI_est_noso_20230820.rda")
pool_noso <- obtain_est(GTest)

# 6) office
load("Main_GI_est_white_20230822.rda")
pool_white <- obtain_est(GTest)

# 7) blue-collar
load("Main_GI_est_blue_20230821.rda")
pool_blue <- obtain_est(GTest)

# 8）fam
load("Main_GI_est_fam_20230820.rda")
pool_fam <- obtain_est(GTest)

pool_res_main <- rbind(pool_fam, pool_carehome, pool_restaurant, pool_noso,
                       pool_blue, pool_white, pool_ICFS, pool_GCLA)

pool_res_main.df <- as.data.frame(pool_res_main)
pool_res_main.df
pool_res_main.df$cluster <- c("Households", "Care homes", "Restaurants", "Nosocomial",
                              "Blue-collar", "Office", "Indoor catering", "General consumption")
pool_res_main.df$cluster <- factor(pool_res_main.df$cluster, levels = 
                                     c("Households", "Care homes", "Restaurants", "Nosocomial",
                                       "Blue-collar", "Office", "Indoor catering", "General consumption"))
pool_res_main.df$X <- seq(1, 8, 1)
colnames(pool_res_main.df)
write.csv(pool_res_main.df, "pooled_GI_est_main.csv")


library(ggplot2)
library(ggpubr)

plotmu <- ggplot() + geom_point(data = pool_res_main.df, aes(x = X, y = mu.pool, color = cluster)) +
  geom_errorbar(data = pool_res_main.df, aes(x = X, ymin = mu.lb, ymax = mu.ub, color = cluster), width = 0.25) +
  scale_color_manual(values = c("#32CD32", "#8B4513", "#FF8C00", "#800080", 
                                "#1E90FF", "#808080", "#FF1493", "#008080" )) +
  scale_x_continuous(breaks = seq(1, 8, 1),
                     labels = pool_res_main.df$cluster) +
  xlab("Cluster setting") + ylab("Days") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "none") +
  ggtitle("Mean generation interval")


plotsd <- ggplot() + geom_point(data = pool_res_main.df, aes(x = X, y = sigma.pool, color = cluster)) +
  geom_errorbar(data = pool_res_main.df, aes(x = X, ymin = sigma.lb, ymax = sigma.ub, color = cluster), width = 0.25) +
  scale_color_manual(values = c("#32CD32", "#8B4513", "#FF8C00", "#800080", 
                                "#1E90FF", "#808080", "#FF1493", "#008080" )) +
  scale_x_continuous(breaks = seq(1, 8, 1),
                     labels = pool_res_main.df$cluster) +
  xlab("Cluster setting") + ylab("Days") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "none") +
  ggtitle("Standard deviation of generation interval")

pool_res_main.df$shape <- pool_res_main.df$mu.pool^2/pool_res_main.df$sigma.pool^2
pool_res_main.df$rate <- pool_res_main.df$mu.pool/pool_res_main.df$sigma.pool^2

plotpdf <- ggplot() + 
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[1],
                                         rate = pool_res_main.df$rate[1])},
                aes(color = paste(1, pool_res_main.df$cluster[1]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[2],
                                         rate = pool_res_main.df$rate[2])},
                aes(color = paste(2, pool_res_main.df$cluster[2]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[3],
                                         rate = pool_res_main.df$rate[3])},
                aes(color = paste(3, pool_res_main.df$cluster[3]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[4],
                                         rate = pool_res_main.df$rate[4])},
                aes(color = paste(4, pool_res_main.df$cluster[4]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[5],
                                         rate = pool_res_main.df$rate[5])},
                aes(color = paste(5, pool_res_main.df$cluster[5]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[6],
                                         rate = pool_res_main.df$rate[6])},
                aes(color = paste(6, pool_res_main.df$cluster[6]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[7],
                                         rate = pool_res_main.df$rate[7])},
                aes(color = paste(7, pool_res_main.df$cluster[7]))) +
  geom_function(fun = function(x){dgamma(x, shape = pool_res_main.df$shape[8],
                                         rate = pool_res_main.df$rate[8])},
                aes(color = paste(8, pool_res_main.df$cluster[8]))) +
  scale_color_manual(
    name = c("Cluster setting"),
    values = c("#32CD32", "#8B4513",
               "#FF8C00", "#800080",
               "#1E90FF", "#808080",
               "#FF1493", "#008080"),
    labels = pool_res_main.df$cluster) +
  scale_x_continuous(limits = c(0, 20)) +
  xlab("Days") + ylab("Probability") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x.bottom = element_line(), 
        axis.line.y.left = element_line(),
        legend.position = c(0.85, 0.5)) +
  ggtitle("Probability density of inferred generation interval distributions")

plotpdf

p1 <- ggarrange(plotmu, plotsd, nrow = 2, ncol = 1, align = "v", labels = c("a", "b"))
pcombine <- ggarrange(p1, plotpdf, nrow = 2, ncol = 1, labels = c("", "c"))
pcombine

ggsave("Fig 3_update.pdf", pcombine, width = 12, height = 9, units = "in")








