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


## infection to report intv
## Gamma mean 7.98 SD 3.3

# shape = mu^2/sd^2 = 7.98^2/3.3^2 = 5.85
# rate = mu/sd^2 = 7.98/3.3^2 = 0.73

infer_Tinfection_alt <- function(clusterdata, intvpara, linelist_ref, iseed){
  set.seed(iseed)
  intvpara <- c(5.85, 0.73)
  allcase <- sort(unique(unlist(clusterdata$cases)))
  ind.allcase <- which(linelist_ref$`HK case no.` %in% allcase)
  reports <- linelist_ref$`Report date`[ind.allcase]
  
  reports.numeric <- as.numeric(reports - as.Date("2020-01-01"))
  sample_infct_to_reprt <- rtrunc(length(ind.allcase),
                                  "gamma",
                                  shape = intvpara[1], rate = intvpara[2],
                                  lb = 1, ub = 21) 
  
  t_infct <- round(reports.numeric - sample_infct_to_reprt)
  
  time_ref <- data.frame(
    case = allcase,
    t_infct = t_infct
  )
  
  infer.infct.clust <- vector("list", length(clusterdata$cluster))
  for(i in 1:length(clusterdata$cluster)){
    case <- unlist(clusterdata$cases[[i]])
    ind.extract <- which(time_ref$case %in% case)
    infer.infct.clust[[i]] <- time_ref$t_infct[ind.extract]
  }
  
  sample_result <- list(
    inferred.Tinfection = infer.infct.clust
  )
  
  return(sample_result)
  
}

sample_times_retail <- vector("list", 100)

progressbar <- txtProgressBar(min = 0, max = 100, style = 3)

tstart <- Sys.time()

for(i in 1:100){
  sample_times_retail[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$retail,
                                               intvpara = c(5.85, 0.73),
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
  sample_times_carehome[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$carehome,
                                                 intvpara = c(5.85, 0.73),
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
  sample_times_office[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$office,
                                               intvpara = c(5.85, 0.73),
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
  sample_times_manual[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$manual,
                                               intvpara = c(5.85, 0.73),
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
  sample_times_noso[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$noso,
                                             intvpara = c(5.85, 0.73),
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
  sample_times_social[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$social,
                                               intvpara = c(5.85, 0.73),
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
  sample_times_restaurant[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$restaurant,
                                                   intvpara = c(5.85, 0.73),
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
  sample_times_household[[i]] <- infer_Tinfection_alt(clusterdata = cluster_n_case$household,
                                                  intvpara = c(5.85, 0.73),
                                                  linelist_ref = linelist_simp,
                                                  iseed = i)
  print(i)
  setTxtProgressBar(progressbar, i)
}
tend <- Sys.time()
tend - tstart


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
source("mixture model.R")
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
sample_times_list_alt <- list(
  household = sample_times_household,
  office = sample_times_office,
  restaurant = sample_times_restaurant,
  manual_labour = sample_times_manual,
  retail = sample_times_retail,
  noso = sample_times_noso,
  social = sample_times_social,
  carehome = sample_times_carehome
)

save(sample_times_list_alt, file = "sample_times_list_alt.rda")
















