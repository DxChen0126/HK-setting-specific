library(mixdist)
library(fitdistrplus)

load("cluster_n_case_list.rda")
load("linelist_ref.rda")

linelist_simp <- linelist_ref

qtrunc <- function(p, spec, lb = -Inf, ub = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(lb, ...) + p*(G(ub, ...) - G(lb, ...)), ...)
  return(tt)
}

rtrunc <- function(n, spec, lb = -Inf, ub = Inf, ...)
{
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, lb = lb, ub = ub,...)
  return(x)
}


# testcluster <- cluster_n_case$household
# allcase <- unlist(testcluster$cases)
# onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase) %>% dplyr::select(`Onset date`)
# reports <- linelist_simp %>% filter(`HK case no.` %in% allcase) %>% dplyr::select(`Report date`)
# 
# onsets.numeric <- as.numeric(onsets$`Onset date` - as.Date("2020-01-01"))
# ind.asymp <- which(is.na(onsets.numeric))
# onsets.numeric <- onsets.numeric[-ind.asymp]
# report.numeric.onset <- as.numeric(reports$`Report date` - as.Date("2020-01-01"))[-ind.asymp]
# t_infect_symp <- estimate_infection_time(onsets.numeric) 
# itr_intv <- report.numeric.onset - t_infect_symp$infection_times
# fitdist(itr_intv, "gamma")

# shape = 6.26 rate = 0.72

####### 
### incu gamma shape 6.25 rate 0.96
### refers to mean 6.5 sd 2.6

infer_Tinfection <- function(clusterdata, incupara, linelist_ref, iseed){
  set.seed(iseed)
  incupara <- c(6.25, 0.96)
  allcase <- sort(unique(unlist(clusterdata$cases)))
  ind.allcase <- which(linelist_ref$`HK case no.` %in% allcase)
  onsets <- linelist_ref$`Onset date`[ind.allcase]
  reports <- linelist_ref$`Report date`[ind.allcase]
  
  onsets.numeric <- as.numeric(onsets - as.Date("2020-01-01"))
  reports.numeric <- as.numeric(reports - as.Date("2020-01-01"))
  
  ind.asymp <- which(is.na(onsets.numeric))
  ind.symp <- which(!is.na(onsets.numeric))
  
  t_infct_symp <- round(onsets.numeric[ind.symp] - 
                          rtrunc(length(ind.symp),
                                 "gamma", shape = incupara[1], rate = incupara[2],
                                 lb = 1, ub = 14))
  
  infct_to_reprt_symp.emp <- reports.numeric[ind.symp] - t_infct_symp
  infct_to_reprt_est <- fitdist(infct_to_reprt_symp.emp, "gamma")
  shape.itv <- infct_to_reprt_est$estimate[1]
  rate.itv <- infct_to_reprt_est$estimate[2]
  itv_boot <- bootdist(infct_to_reprt_est, niter = 100)
  mu.itv <- shape.itv/rate.itv
  sd.itv <- sqrt(shape.itv/rate.itv^2)
  mu.itv.boot <- itv_boot$estim[,1]/itv_boot$estim[,2]
  sd.itv.boot <- sqrt(itv_boot$estim[,1]/itv_boot$estim[,2]^2)
  shape.itv.sd <- infct_to_reprt_est$sd[1]
  rate.itv.sd <- infct_to_reprt_est$sd[2]
  
  sample_infct_to_reprt_asymp <- rtrunc(length(ind.asymp),
                                        "gamma",
                                        shape = shape.itv, rate = rate.itv,
                                        lb = 1, ub = 21)
  
  t_infct_asymp <- round(reports.numeric[ind.asymp] - sample_infct_to_reprt_asymp)
  time_ref <- data.frame(
    case = c(allcase[ind.symp], allcase[ind.asymp]),
    t_infct = c(t_infct_symp, t_infct_asymp)
  )
  
  infer.infct.clust <- vector("list", length(clusterdata$cluster))
  for(i in 1:length(clusterdata$cluster)){
    case <- unlist(clusterdata$cases[[i]])
    ind.extract <- which(time_ref$case %in% case)
    infer.infct.clust[[i]] <- time_ref$t_infct[ind.extract]
  }
  
  sample_result <- list(
    intv_est = list(
      mu = as.numeric(mu.itv),
      sd = as.numeric(sd.itv),
      mu.boot = mu.itv.boot,
      sd.boot = sd.itv.boot
    ),
    inferred.Tinfection = infer.infct.clust
  )
  
  return(sample_result)

}


infer_Tinfection(clusterdata = cluster_n_case$restaurant,
                 incupara = c(6.25, 0.96),
                 linelist_ref = linelist_ref,
                 iseed = 1)

cluster_n_case$restaurant$cases[[147]]

###### it works

sample_times_retail <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_retail[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$retail,
                                               incupara = c(6.25, 0.96),
                                               linelist_ref = linelist_simp,
                                               iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

sample_times_carehome <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_carehome[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$carehome,
                                               incupara = c(6.25, 0.96),
                                               linelist_ref = linelist_simp,
                                               iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

sample_times_office <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_office[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$office,
                                                 incupara = c(6.25, 0.96),
                                                 linelist_ref = linelist_simp,
                                                 iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

sample_times_manual <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_manual[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$manual,
                                               incupara = c(6.25, 0.96),
                                               linelist_ref = linelist_simp,
                                               iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

sample_times_noso <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_noso[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$noso,
                                               incupara = c(6.25, 0.96),
                                               linelist_ref = linelist_simp,
                                               iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart


sample_times_social <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_social[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$social,
                                             incupara = c(6.25, 0.96),
                                             linelist_ref = linelist_simp,
                                             iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

sample_times_restaurant <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_restaurant[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$restaurant,
                                               incupara = c(6.25, 0.96),
                                               linelist_ref = linelist_simp,
                                               iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

sample_times_household <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_household[[i]] <- infer_Tinfection(clusterdata = cluster_n_case$household,
                                                   incupara = c(6.25, 0.96),
                                                   linelist_ref = linelist_simp,
                                                   iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

delay_intv_est <- function(inputdata){
  mu.pool = mean(unlist(lapply(inputdata, function(x) x$intv_est$mu)))
  sd.pool = mean(unlist(lapply(inputdata, function(x) x$intv_est$sd)))
  mu.boot.pool = unlist(lapply(inputdata, function(x) x$intv_est$mu.boot))
  sd.boot.pool = unlist(lapply(inputdata, function(x) x$intv_est$sd.boot))
  return(data.frame(
    mu_pool = mu.pool,
    mu_LB = quantile(mu.boot.pool, 0.025),
    mu_UB = quantile(mu.boot.pool, 0.975),
    sd_pool = sd.pool,
    sd_LB = quantile(sd.boot.pool, 0.025),
    sd_UB = quantile(sd.boot.pool, 0.975)
  ))
}

delay_intv_household <- delay_intv_est(sample_times_list$household)
delay_intv_office <- delay_intv_est(sample_times_list$office)
delay_intv_restaurant <- delay_intv_est(sample_times_list$restaurant)
delay_intv_manual <- delay_intv_est(sample_times_list$manual_labour)
delay_intv_retail <- delay_intv_est(sample_times_list$retail)
delay_intv_noso <- delay_intv_est(sample_times_list$noso)
delay_intv_social <- delay_intv_est(sample_times_list$social)
delay_intv_carehome <- delay_intv_est(sample_times_list$carehome)

delay_intv_list <- list(
  household = delay_intv_household,
  office = delay_intv_office,
  restaurant = delay_intv_restaurant,
  manual_labour = delay_intv_manual,
  retail = delay_intv_retail,
  noso = delay_intv_noso,
  social = delay_intv_social,
  carehome = delay_intv_carehome
)

#save(delay_intv_list, file = "infct_reprt_delay_list.rda")
load("infct_reprt_delay_list.rda")
delay_intv_list


arrange_ICCI <- function(data){
  n <- length(data)
  temp_list <- vector("list", n)
  for(i in 1:n){
    temp_data_sort <- sort(data[[i]])
    if(length(temp_data_sort) >= 12){
      temp_data_sort <- sort(sample(temp_data_sort, 11, replace = FALSE))
    }
    ICCI <- (temp_data_sort - temp_data_sort[1])[-1]
    temp_list[[i]] <- ICCI
  }
  return(unlist(temp_list))
}

######### adjust sampled times

### adjust office cluster
office.subtract <- read_excel("subtract_office_case.xlsx")
ind.check <- which(cluster_n_case$office$cluster %in% office.subtract$cluster_name)
ind.case.subtract <- vector("list", 4)
for(i in 1:4){
  ind.case.subtract[[i]] <- which(cluster_n_case$office$cases[[ind.check[i]]] 
                                  %in% office.subtract$subtract_case[i])
}

for(i in 1:100){
  for(j in 1:4){
    sample_times_office[[i]]$inferred.Tinfection[[ind.check[j]]] <-
      sample_times_office[[i]]$inferred.Tinfection[[ind.check[j]]][-unlist(ind.case.subtract[j])]
  }
}

testICCI <- arrange_ICCI(sample_times_office[[2]]$inferred.Tinfection)

GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)

retail.subtract <- read_excel("subtract_retail_case.xlsx")
ind.check <- which(cluster_n_case$retail$cluster %in% retail.subtract$cluster_name)
ind.case.subtract  <- which(cluster_n_case$retail$cases[[ind.check]] 
                                  %in% retail.subtract$subtract_case)

for(i in 1:100){
    sample_times_retail[[i]]$inferred.Tinfection[ind.check] <-
      sample_times_retail[[i]]$inferred.Tinfection[ind.check][-ind.case.subtract]
}

testICCI <- arrange_ICCI(sample_times_retail[[100]]$inferred.Tinfection)

GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)

noso.subtract <- read_excel("subtract_noso_case.xlsx")
ind.check <- which(cluster_n_case$noso$cluster %in% noso.subtract$cluster_name)
ind.case.subtract  <- which(cluster_n_case$noso$cases[[ind.check]] 
                            %in% noso.subtract$subtract_case)

for(i in 1:100){
  sample_times_noso[[i]]$inferred.Tinfection[ind.check] <-
    sample_times_noso[[i]]$inferred.Tinfection[ind.check][-ind.case.subtract]
}

testICCI <- arrange_ICCI(sample_times_noso[[1]]$inferred.Tinfection)

GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)

manual.subtract <- read_excel("subtract_manual_case.xlsx")
ind.check <- which(cluster_n_case$manual$cluster %in% manual.subtract$cluster_name)
ind.case.subtract <- vector("list", 3)
for(i in 1:3){
  ind.case.subtract[[i]] <- which(cluster_n_case$manual$cases[[ind.check[i]]] 
                                  %in% manual.subtract$subtract_case[i])
}

for(i in 1:100){
  for(j in 1:3){
    sample_times_manual[[i]]$inferred.Tinfection[[ind.check[j]]] <-
      sample_times_manual[[i]]$inferred.Tinfection[[ind.check[j]]][-unlist(ind.case.subtract[j])]
  }
}

testICCI <- arrange_ICCI(sample_times_manual[[10]]$inferred.Tinfection)

GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)

social.subtract <- read_excel("subtract_social_case.xlsx")
ind.check <- which(cluster_n_case$social$cluster %in% social.subtract$cluster_name)
ind.case.subtract  <- which(cluster_n_case$social$cases[[ind.check]] 
                            %in% as.numeric(unlist(str_split(social.subtract$subtract_case, pattern = "; "))))

for(i in 1:100){
  sample_times_social[[i]]$inferred.Tinfection[ind.check] <-
    sample_times_social[[i]]$inferred.Tinfection[ind.check][-ind.case.subtract]
}

testICCI <- arrange_ICCI(sample_times_social[[10]]$inferred.Tinfection)
GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)

restaurant.subtract <- read_excel("subtract_restaurant_case.xlsx")
ind.check <- which(cluster_n_case$restaurant$cluster %in% restaurant.subtract$cluster_name)
ind.case.subtract <- vector("list", 9)
for(i in 1:9){
  ind.case.subtract[[i]] <- which(cluster_n_case$restaurant$cases[[ind.check[i]]] 
                                  %in% as.numeric(unlist(str_split(restaurant.subtract$subtract_case, pattern = "; "))))
}

for(i in 1:100){
  for(j in 1:9){
    sample_times_restaurant[[i]]$inferred.Tinfection[[ind.check[j]]] <-
      sample_times_restaurant[[i]]$inferred.Tinfection[[ind.check[j]]][-unlist(ind.case.subtract[j])]
  }
}

testICCI <- arrange_ICCI(sample_times_restaurant[[10]]$inferred.Tinfection)
GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)

household.subtract <- read_excel("subtract_household_case.xlsx")
ind.check <- which(cluster_n_case$household$cluster %in% household.subtract$cluster_name)
ind.case.subtract <- vector("list", 34)
for(i in 1:34){
  ind.case.subtract[[i]] <- which(cluster_n_case$household$cases[[ind.check[i]]] 
                                  %in% as.numeric(unlist(str_split(household.subtract$subtract_case, pattern = "; "))))
}

for(i in 1:100){
  for(j in 1:34){
    sample_times_household[[i]]$inferred.Tinfection[[ind.check[j]]] <-
      sample_times_household[[i]]$inferred.Tinfection[[ind.check[j]]][-unlist(ind.case.subtract[j])]
  }
}

testICCI <- arrange_ICCI(sample_times_household[[10]]$inferred.Tinfection)
GI_mix_est(data = testICCI, N = 50, startmu = 5, startsig = 5)


############# save these 
sample_times_list <- list(
  household = sample_times_household,
  office = sample_times_office,
  restaurant = sample_times_restaurant,
  manual_labour = sample_times_manual,
  retail = sample_times_retail,
  noso = sample_times_noso,
  social = sample_times_social,
  carehome = sample_times_carehome
)

save(sample_times_list, file = "sample_times_list.rda")

#######
########## given another descriptive of asymp proportion in each cluster

load("cluster_n_case_list.rda")
load("linelist_ref.rda")
load("clustersize_list.rda")

allcase_household <- as.numeric(c(unlist(cluster_n_case$household$cases),
    clustersize_list$household$clustername[clustersize_list$household$clustersize == 1]))
household_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_household) %>% 
  dplyr::select(`Onset date`)
asymp_prop_household <- length(which(is.na(household_onsets)))/nrow(household_onsets)
# 0.2093535

allcase_office <- as.numeric(c(unlist(cluster_n_case$office$cases),
                                  clustersize_list$office$clustername[clustersize_list$office$clustersize == 1]))
office_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_office) %>% 
  dplyr::select(`Onset date`)
asymp_prop_office <- length(which(is.na(office_onsets)))/nrow(office_onsets)
# 0.1445312

allcase_restaurant <- as.numeric(c(unlist(cluster_n_case$restaurant$cases),
                               clustersize_list$restaurant$clustername[clustersize_list$restaurant$clustersize == 1]))
allcase_restaurant <- allcase_restaurant[which(!is.na(allcase_restaurant))]
restaurant_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_restaurant) %>% 
  dplyr::select(`Onset date`)
asymp_prop_restaurant <- length(which(is.na(restaurant_onsets)))/nrow(restaurant_onsets)
# 0.1749664

allcase_manual <- as.numeric(c(unlist(cluster_n_case$manual$cases),
                                   clustersize_list$manual$clustername[clustersize_list$manual$clustersize == 1]))
manual_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_manual) %>% 
  dplyr::select(`Onset date`)
asymp_prop_manual <- length(which(is.na(manual_onsets)))/nrow(manual_onsets)
# 0.2627866

allcase_retail <- as.numeric(c(unlist(cluster_n_case$retail$cases),
                               clustersize_list$retail$clustername[clustersize_list$retail$clustersize == 1]))
retail_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_retail) %>% 
  dplyr::select(`Onset date`)
asymp_prop_retail <- length(which(is.na(retail_onsets)))/nrow(retail_onsets)
# 0.1151832

allcase_noso <- as.numeric(c(unlist(cluster_n_case$noso$cases),
                               clustersize_list$noso$clustername[clustersize_list$noso$clustersize == 1]))
noso_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_noso) %>% 
  dplyr::select(`Onset date`)
asymp_prop_noso <- length(which(is.na(noso_onsets)))/nrow(noso_onsets)
# 0.1550802

allcase_social <- as.numeric(c(unlist(cluster_n_case$social$cases),
                             clustersize_list$social$clustername[clustersize_list$social$clustersize == 1]))
social_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_social) %>% 
  dplyr::select(`Onset date`)
asymp_prop_social <- length(which(is.na(social_onsets)))/nrow(social_onsets)
# 0.2816901

allcase_care <- as.numeric(c(unlist(cluster_n_case$carehome$cases),
                               clustersize_list$carehome$clustername[clustersize_list$care$clustersize == 1]))
care_onsets <- linelist_simp %>% filter(`HK case no.` %in% allcase_care) %>% 
  dplyr::select(`Onset date`)
asymp_prop_care <- length(which(is.na(care_onsets)))/nrow(care_onsets)
# 0.39375

########### And then apply the GI mixture model for estimation
source("mixture model.R")

load("sample_times_list.rda")
library(doParallel)

ncores <- detectCores()


testICCI <- arrange_ICCI(sample_times_restaurant[[1]]$inferred.Tinfection)

GI_mix_est(data = testICCI, N = 40, startmu = 10, startsig = 10)

save(allclustcomb.df, file = "all_cluster_comb_df.rda")


sample_times_noso <- sample_times_list$noso

GI_est_noso <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100){
  set.seed(i)
  ICCI <- arrange_ICCI(sample_times_noso[[i]]$inferred.Tinfection)
  ICCI_est <- GI_mix_est(data = ICCI, N = 50, startmu = 5, startsig = 5)
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  estboot <- foreach(i = 1:100, .combine = 'append') %dopar% {
    GI_mix_est(data = sample(ICCI, length(ICCI), replace = T),
               N = 50, startmu = 5, startsig = 5)
  }
  stopCluster(myCluster)
  GI_est_noso[[i]] <- list(est = ICCI_est,
                             est_boot = estboot)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

save(GI_est_noso, file = "GI_est_noso.rda")

load("GI_est_noso.rda")


lapply(GI_est_noso, function(x) x$est[1]) %>% unlist() %>% mean()
# 6.885313
lapply(GI_est_noso, function(x) x$est_boot[seq(1, 500, 5)]) %>% unlist() %>% quantile(c(0.025, 0.975))
# 3.853798 11.737301 

lapply(GI_est_noso, function(x) x$est[2]) %>% unlist() %>% mean()
# 2.523112
lapply(GI_est_noso, function(x) x$est_boot[seq(2, 500, 5)]) %>% unlist() %>% quantile(c(0.025, 0.975))
# 1.010976 5.003381  

rm(GI_est_noso)

sample_times_care <- sample_times_list$carehome

GI_est_care <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100){
  set.seed(i)
  ICCI <- arrange_ICCI(sample_times_care[[i]]$inferred.Tinfection)
  ICCI_est <- GI_mix_est(data = ICCI, N = 50, startmu = 5, startsig = 5)
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  estboot <- foreach(i = 1:100, .combine = 'append') %dopar% {
    GI_mix_est(data = sample(ICCI, length(ICCI), replace = T),
               N = 50, startmu = 5, startsig = 5)
  }
  stopCluster(myCluster)
  GI_est_care[[i]] <- list(est = ICCI_est,
                           est_boot = estboot)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

save(GI_est_care, file = "GI_est_care.rda")
load("GI_est_care.rda")

lapply(GI_est_care, function(x) x$est[1]) %>% unlist() %>% mean()
#  7.025881
lapply(GI_est_care, function(x) x$est_boot[seq(1, 500, 5)]) %>% unlist() %>% quantile(c(0.025, 0.975))
# 3.73775 10.32031 

lapply(GI_est_care, function(x) x$est[2]) %>% unlist() %>% mean()
# 2.062089
lapply(GI_est_care, function(x) x$est_boot[seq(2, 500, 5)]) %>% unlist() %>% quantile(c(0.025, 0.975))
# 1.088108 4.349725 

################
sample_times_restaurant <- sample_times_list$restaurant

GI_est_restaurant <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100){
  set.seed(i)
  ICCI <- arrange_ICCI(sample_times_restaurant[[i]]$inferred.Tinfection)
  ICCI_est <- GI_mix_est(data = ICCI, N = 50, startmu = 5, startsig = 5)
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  estboot <- foreach(i = 1:100, .combine = 'append') %dopar% {
    GI_mix_est(data = sample(ICCI, length(ICCI), replace = T),
               N = 50, startmu = 5, startsig = 5)
  }
  stopCluster(myCluster)
  GI_est_restaurant[[i]] <- list(est = ICCI_est,
                           est_boot = estboot)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart

save(GI_est_restaurant, file = "GI_est_restaurant.rda")
load("GI_est_restaurant.rda")

lapply(GI_est_restaurant, function(x) x$est[1]) %>% unlist() %>% mean()
#  7.025881
lapply(GI_est_restaurant, function(x) x$est_boot[seq(1, 500, 5)]) %>% unlist() %>% quantile(c(0.025, 0.975))
# 3.73775 10.32031 

lapply(GI_est_restaurant, function(x) x$est[2]) %>% unlist() %>% mean()
# 2.062089
lapply(GI_est_restaurant, function(x) x$est_boot[seq(2, 500, 5)]) %>% unlist() %>% quantile(c(0.025, 0.975))
# 1.088108 4.349725 

### and so for other settings...





