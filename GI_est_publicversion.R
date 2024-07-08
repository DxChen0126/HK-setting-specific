rm(list = ls())

library(dplyr)
library(foreach)
library(doParallel)

ncores <- detectCores()

load("main_GI_100times_sample_20230803.rda")

source("mixture model 0801.R")


progressbar <- txtProgressBar(min = 0, max = 100, style = 3)
tstart <- Sys.time()
SampleGT <- sampleGT_record[[8]]     ### index 1 to 8 correspond to each cluster setting
GTest <- vector("list", 100)
for(j in 1:100){
  set.seed(j)
  tempdata <- SampleGT[[1 + (j-1)*2]]
  tempest <- GI_mix_est_revise(data = tempdata, N = 100, startmu = 10, startsig = 10)
  # save(tempest, file = paste0("./sample infection record/GCLA/GI_est_GCLA_", j, "th_est.rda"))
  # do parallel for the bootstap part
  myCluster <- makeCluster(ncores - 1, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  estboot <- foreach(i = 1:100, .combine = 'append',
                     .packages=c("dplyr"), .errorhandling = 'pass') %dopar% {
                       GI_mix_est_revise(data = sample(tempdata, length(tempdata), replace = T),
                                         N = 100, startmu = 10, startsig = 10)
                     }
  # save(estboot, file = paste0("./sample infection record/GCLA/GI_boot_GCLA_", j, "th_sample.rda"))
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

### and so for each cluster setting


### a function to obtain estimation result

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

### then do for each cluster setting

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






