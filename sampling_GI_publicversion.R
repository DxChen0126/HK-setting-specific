### 1) for each cluster type's symptomatic case, sample from 
# onset to infection by IP dist (fixed) then fit dist from infection to
# report, repeat sampling for 100 times to provide dist from infection to report

rm(list = ls())
load("para_dist.rda")
load("cluster_related_cases.rda")
load("cluster_uniqID.rda")
load("linelist_ref.rda")
load("clustIDlist_ref.rda")
library(fitdistrplus)
library(dplyr)

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

### the goal is to fit dist of infection to report intv based on symptomatic cases in each cluster
dist_par_list
incu.shape <- dist_par_list$Incu_gamma[1]
incu.rate <- dist_par_list$Incu_gamma[2]

fit_inft_rpt_dist <- function(datainput, incushape, incurate, iseed){
  # datainput is the clustdata (contains cluster names)
  # to make the results reproduceable
  set.seed(iseed)
  nclust <- length(datainput)
  inft_rpt_intv_list <- vector("list", nclust)
  clustdata = datainput
  for(j in 1:nclust){
    refID <- clustdata[j]
    ind.ID <- which(clustall.uniq == refID)
    refcase <- clustall.uniq.relatecase.alt[[ind.ID]]
    datatmp <- linelist_simp %>% filter(`HK case no.` %in% refcase)
    ind.symp <- which(!is.na(datatmp$`Onset date`))
    incu.intv <- rtrunc(length(ind.symp), "gamma", lb = 1, ub = 20,
                        shape = incushape, rate = incurate)
    inft_rpt_intv <- (datatmp[ind.symp,]$`Report date` - datatmp[ind.symp,]$`Onset date`) %>% 
      as.numeric() + incu.intv %>% round()
    # in case of negative intervals
    inft_rpt_intv[inft_rpt_intv <= 0] <- 1
    inft_rpt_intv_list[[j]] <- inft_rpt_intv
  }
  intv <- unlist(inft_rpt_intv_list)
  myfit <- fitdist(intv, "gamma", method = "mle")
  return(list(shape.est = myfit$estimate[1],
              rate.est = myfit$estimate[2],
              shape.se = myfit$sd[1],
              rate.se = myfit$sd[2]))
}

library(foreach)
library(doParallel)

ncores <- detectCores()

dist_fit <- vector("list", 8)

progressbar <- txtProgressBar(min = 0, max = 8, style = 3)

tstart <- Sys.time()
for(m in 1:8){
  clustdata = clustID_list[[m]]
  # set 100 sampling times
  nsample = 100
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  dist_fit[[m]] <- foreach(i = 1:nsample, .combine = 'append',
                           .packages=c("dplyr", "fitdistrplus")) %dopar% {
                             fit_inft_rpt_dist(datainput = clustdata,
                                               incushape = incu.shape,
                                               incurate = incu.rate,
                                               iseed = i)
                           }
  stopCluster(myCluster)
  print(m)
  setTxtProgressBar(progressbar, m)
}

tend <- Sys.time()
tend - tstart

inft_rpt_fitdist <- vector("list", 8)
for(i in 1:8){
  restmp <- dist_fit[[i]]
  shape.estpool <- restmp[seq(1, 400, 4)] %>% unlist() %>% mean()
  rate.estpool <- restmp[seq(2, 400, 4)] %>% unlist() %>% mean()
  shape.var <- restmp[seq(1, 400, 4)] %>% unlist() %>% var() + mean(unlist(restmp[seq(3, 400, 4)])^2)
  rate.var <- restmp[seq(2, 400, 4)] %>% unlist() %>% var() + mean(unlist(restmp[seq(4, 400, 4)])^2)
  para.res <- c(shape.estpool, shape.estpool - 1.96 * sqrt(shape.var), 
                shape.estpool + 1.96 * sqrt(shape.var),
                rate.estpool, rate.estpool - 1.96 * sqrt(rate.var),
                rate.estpool + 1.96 * sqrt(rate.var))
  names(para.res) <- c("shape_poolest", "shape_LB", "shape_UB",
                       "rate_poolest", "rate_LB", "rate_UB")
  inft_rpt_fitdist[[i]] <- para.res
}

inft_rpt_mu <- inft_rpt_sd <- numeric(8)
for(i in 1:8){
  inft_rpt_mu[i] <- inft_rpt_fitdist[[i]][1]/inft_rpt_fitdist[[i]][4]
  inft_rpt_sd[i] <- sqrt(inft_rpt_fitdist[[i]][1]/inft_rpt_fitdist[[i]][4]^2)
}
inft_rpt_mu
# 11.080059  9.139358 11.741738 11.055049 10.933932 11.169689 10.290082 12.733400
inft_rpt_sd
# 3.994394 3.278436 4.083453 4.019594 4.025372 3.737219 3.617759 3.906562

# compare with that NC paper:
dist_par_list[[2]][1]/dist_par_list[[2]][2] # 7.979609 
sqrt(dist_par_list[[2]][1]/dist_par_list[[2]][2]^2) # 3.302913 

save(inft_rpt_fitdist, file = "infection_to_report_dist_bycluster.rda")

## now we try to sample the infection date for each cluster's asymptomatic cases
# based on such distribution stratified by cluster types
#function for inferring infection times
infection_time_infer <- function(inputdata, incu_shape, incu_rate,
                                 Ditt_shape, Ditt_rate){
  ind.asymp <- which(is.na(inputdata$`Onset date`))
  ind.symp <- which(!is.na(inputdata$`Onset date`))
  infect_time <- c()
  # 1) no asymptomatic cases
  if(length(ind.asymp)==0){
    nsymp <- nrow(inputdata)
    intv.symp <- rtrunc(nsymp, "gamma", lb = 1, ub = 20,
                        shape = incu_shape, rate = incu_rate)
    infect_time <- inputdata$`Onset date` - intv.symp
  }
  # 2) no symptomatic cases
  if(length(ind.symp)==0){
    nasymp <- nrow(inputdata)
    intv.asymp <- rtrunc(nasymp, "gamma", lb = 2, ub = 20,
                         shape = Ditt_shape, rate = Ditt_rate)
    infect_time <- inputdata$`Report date` - intv.asymp
  }
  # 3) some asymptomatic some symptomatic
  if(length(ind.symp)>0 & length(ind.asymp)>0){
    nsymp <- length(ind.symp)
    nasymp <- length(ind.asymp)
    intv.symp <- rtrunc(nsymp, "gamma", lb = 1, ub = 20,
                        shape = incu_shape, rate = incu_rate)
    intv.asymp <- rtrunc(nasymp, "gamma", lb = 2, ub = 20,
                         shape = Ditt_shape, rate = Ditt_rate)
    infect_time <- c((inputdata[ind.asymp,]$`Report date` - intv.asymp),
                     (inputdata[ind.symp,]$`Onset date` - intv.symp))
  }
  
  return(infect_time)
}

sample_infect_intv <- function(datainput, incu_shape, incu_rate,
                               Ditt_shape, Ditt_rate, iseed){
  # input data is a vector that contains all clusters in each cluster category
  #i.e. clust_famID <- clustID_list$fam
  set.seed(iseed)
  # incu_shape = dist_par_list$Incu_gamma[1]
  # incu_rate = dist_par_list$Incu_gamma[2]
  # Ditt_shape = dist_par_list$Ditt_gamma[1]
  # Ditt_rate = dist_par_list$Ditt_gamma[2]
  # datainput <- clustdata
  nclust <- length(datainput)
  infect_intv_list <- vector("list", nclust)
  for(j in 1:nclust){
    refID <- datainput[j]
    ind.ID <- which(clustall.uniq == refID)
    refcase <- clustall.uniq.relatecase.alt[[ind.ID]]
    datatmp <- linelist_simp %>% filter(`HK case no.` %in% refcase)
    infect_times <- infection_time_infer(datatmp, 
                                         incu_shape = incu_shape, incu_rate = incu_rate,
                                         Ditt_shape = Ditt_shape, Ditt_rate = Ditt_rate)
    # prepare for the mixture model, cannot deal with large size
    if(length(infect_times) > 11){
      infect_times <- sample(infect_times, 11, replace = FALSE)
    }
    # create infection intervals
    infect_times <- sort(infect_times) # sort from earliest to latest
    infect_intv <- round(as.numeric(infect_times[-1] - infect_times[1])) 
    infect_intv_list[[j]] <- infect_intv
  }
  
  allinfect_intv <- unlist(infect_intv_list)
  # further drop out infection interval > 30, which looks highly impossible
  allinfect_intv <- allinfect_intv[allinfect_intv <= 30]
  return(list(
    allintervals = allinfect_intv,
    totalsize = length(allinfect_intv)
  ))
}

# for each cluster type, sample infection intervals for 100 times, keep the sampling
# record, then conduct the estimation later
library(foreach)
library(doParallel)

ncores <- detectCores()

sampleGT_record <- vector("list", 8)

progressbar <- txtProgressBar(min = 0, max = 8, style = 3)
tstart <- Sys.time()
for(m in 1:8){
  clustdata = clustID_list[[m]]
  incu.shape = dist_par_list$Incu_gamma[1]
  incu.rate = dist_par_list$Incu_gamma[2]
  inft_rpt.shape = inft_rpt_fitdist[[m]][1]
  inft_rpt.rate = inft_rpt_fitdist[[m]][4]
  # set 100 sampling times
  nsample = 100
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  sampleGT_record[[m]] <- foreach(i = 1:nsample, .combine = 'append',
                                  .packages=c("dplyr")) %dopar% {
                                    sample_infect_intv(datainput = clustdata,
                                                       incu_shape = incu.shape,
                                                       incu_rate = incu.rate,
                                                       Ditt_shape = inft_rpt.shape,
                                                       Ditt_rate = inft_rpt.rate,
                                                       iseed = i)
                                  }
  stopCluster(myCluster)
  print(m)
  setTxtProgressBar(progressbar, m)
}

tend <- Sys.time()
tend - tstart
save(sampleGT_record, file = "./sample infection record/main_GI_100times_sample_20230803.rda")
