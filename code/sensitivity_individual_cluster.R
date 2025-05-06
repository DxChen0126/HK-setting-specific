### these two are the anonymized data contained cluster name, cases involved and 
### cluster size
load("cluster_n_case_list.rda")
load("clustersize_list.rda")

### this is the processed dataframe for visualization and model estimation
load("allclustcomb_df.rda")

library(tidyverse)
library(readxl)
library(openxlsx)

load("linelist_ref.rda")
# build a function to define clusters as index (first onset) + 8 days
# asymptomatic use report date
# to determine cluster size

linelist_symp <- linelist_ref %>% 
  mutate(onset_alt = case_when(
    is.na(`Onset date`) ~ `Report date`,
    !is.na(`Onset date`) ~ `Onset date`
  ))

View(linelist_symp)

alt_size <- function(case_ids, linelist_symp, SIthres){
#  case_ids <- cluster_n_case$restaurant$cases[[10]]
  n_case <- length(case_ids)
  onsets <- linelist_symp %>% filter(`HK case no.` %in% case_ids) %>%
    dplyr::select(onset_alt)
  onset_order <- sort(onsets$onset_alt)
  onset_length <- onset_order[n_case] -  onset_order[1]
  if(onset_length <= SIthres){
    n_size <- n_case
  }
  if(onset_length > SIthres){
    # make subset of cluster sizes
    n_sub <- ceiling((onset_length + 1) / (SIthres + 1))
    n_size <- numeric(n_sub)
    for(i in 1:n_sub){
      n_size[i] <- length(which(onset_order >= onset_order[1] + (i - 1) * (SIthres + 1) &
                    onset_order <= onset_order[1] + (i - 1) * (SIthres + 1) + SIthres))
    }
  }
  return(n_size)
}

# test function
alt_size(cluster_n_case$household$cases[[3]], 
         linelist_symp, SIthres = 5)
# works

# check household cluster
length(cluster_n_case$household$cluster) # 1782

household_size_alt <- vector("list", 1782)
for(j in 1:1782){
  household_size_alt[[j]] <- alt_size(cluster_n_case$household$cases[[j]],
                                      linelist_symp, SIthres = 7)
  print(j)
}

View(household_size_alt)

cluster_n_case$household$size_alt <- household_size_alt

####### check the rest of the other clusters
length(cluster_n_case$restaurant$cases)

restaurant_size_alt <- vector("list", 175)
for(j in 1:175){
  restaurant_size_alt[[j]] <- alt_size(cluster_n_case$restaurant$cases[[j]],
                                      linelist_symp, SIthres = 7)
  print(j)
}

View(restaurant_size_alt)

cluster_n_case$restaurant$size_alt <- restaurant_size_alt

length(cluster_n_case$office$cases)

office_size_alt <- vector("list", 88)
for(j in 1:88){
  office_size_alt[[j]] <- alt_size(cluster_n_case$office$cases[[j]],
                                       linelist_symp, SIthres = 7)
  print(j)
}

View(office_size_alt)

cluster_n_case$office$size_alt <- office_size_alt

length(cluster_n_case$manual$cases)

manual_size_alt <- vector("list", 75)
for(j in 1:75){
  manual_size_alt[[j]] <- alt_size(cluster_n_case$manual$cases[[j]],
                                   linelist_symp, SIthres = 8)
  print(j)
}

View(manual_size_alt)

cluster_n_case$manual$size_alt <- manual_size_alt

length(cluster_n_case$noso$cases)

noso_size_alt <- vector("list", 34)
for(j in 1:34){
  noso_size_alt[[j]] <- alt_size(cluster_n_case$noso$cases[[j]],
                                   linelist_symp, SIthres = 9)
  print(j)
}

View(noso_size_alt)

cluster_n_case$noso$size_alt <- noso_size_alt

length(cluster_n_case$retail$cases)

retail_size_alt <- vector("list", 17)
for(j in 1:17){
  retail_size_alt[[j]] <- alt_size(cluster_n_case$retail$cases[[j]],
                                 linelist_symp, SIthres = 8)
  print(j)
}

View(retail_size_alt)

cluster_n_case$retail$size_alt <- retail_size_alt



length(cluster_n_case$social$cases) # 34

social_size_alt <- vector("list", 34)
for(j in 1:34){
  social_size_alt[[j]] <- alt_size(cluster_n_case$social$cases[[j]],
                                   linelist_symp, SIthres = 8)
  print(j)
}

View(social_size_alt)
cluster_n_case$social$size_alt <- social_size_alt

length(cluster_n_case$carehome$cases) # 23

carehome_size_alt <- vector("list", 23)
for(j in 1:23){
  carehome_size_alt[[j]] <- alt_size(cluster_n_case$carehome$cases[[j]],
                                   linelist_symp, SIthres = 9)
  print(j)
}

View(carehome_size_alt)
cluster_n_case$carehome$size_alt <- carehome_size_alt

save(cluster_n_case, file = "cluster_n_case_20250227.rda")

source("tupper.R")

######### household setting

householdsize <- unlist(household_size_alt)
householdsize <- householdsize[householdsize != 0]

householdsize <- c(householdsize, rep(1, 1536))

# individual ascertainment, q2 = 1, give q1 a value
household.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.9, 
                                                   q2 = 1, 
                                                   clusters = householdsize))

outputres <- function(modelres){
  data.frame(
    dispersion_k = modelres$a_n_b.est[1],
    meanclust = 1 + modelres$a_n_b.est[1] * modelres$a_n_b.est[2],
    Rc = modelres$a_n_b.est[1] * modelres$a_n_b.est[2]
  )
}

outputres(household.mod1)
# dispersion_k meanclust       Rc
# 0.973806    1.7671 0.7670996
conf.household.mod1 <- conf_ellipse(inputdata = list(a.est = household.mod1$a_n_b.est[1],
                                                    b.est = household.mod1$a_n_b.est[2],
                                                    hessmat = household.mod1$hessian_mat),
                                   conflevel = 0.95,
                                   npoint = 100)
range(conf.household.mod1$x.ellipsepoints) # 1.720501 1.813698
range(conf.household.mod1$y.ellipsepoints) # 0.8063177 1.1412943


household.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 1, 
                                                   q2 = 0.9, 
                                                   clusters = householdsize))


outputres(household.mod2)
# dispersion_k meanclust       Rc
# 0.8618233  1.648818 0.6488181
conf.household.mod2 <- conf_ellipse(inputdata = list(a.est = household.mod2$a_n_b.est[1],
                                                     b.est = household.mod2$a_n_b.est[2],
                                                     hessmat = household.mod2$hessian_mat),
                                    conflevel = 0.95,
                                    npoint = 100)
range(conf.household.mod2$x.ellipsepoints) # 1.608484 1.689152
range(conf.household.mod2$y.ellipsepoints) # 0.7224545 1.0011922

household.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.9, 
                                                   q2 = 0.9, 
                                                   clusters = householdsize))

outputres(household.mod3)
# dispersion_k meanclust       Rc
# 0.8590911  1.718889 0.7188888
conf.household.mod3 <- conf_ellipse(inputdata = list(a.est = household.mod3$a_n_b.est[1],
                                                     b.est = household.mod3$a_n_b.est[2],
                                                     hessmat = household.mod3$hessian_mat),
                                    conflevel = 0.95,
                                    npoint = 100)
range(conf.household.mod3$x.ellipsepoints) # 1.674175 1.763602
range(conf.household.mod3$y.ellipsepoints) # 0.7167691 1.0014130



##########COLABORATE WITH Daniel's suggestion to use the update sampling
########## of w1 w2 w3 to have an idea of each generation's cluster size


load("sample_times_list.rda")
View(sample_times_list$household[[1]]$inferred.Tinfection)

load("GI_est_household.rda")
View(GI_est_household[[1]]$est)

source("mixture model.R")

check_likelihood <- function(x, GIest){
  mu = GIest[1]
  sigma = GIest[2]
  w1 = GIest[3]
  w2 = GIest[4] 
  w3 = GIest[5] 
  w4 = 1 - w1 - w2 - w3
  
  a.pri = (mu^2)/(sigma^2) 
  b.pri = mu/(sigma^2)
  
  a.sec = ((2*mu)^2)/(2*sigma^2) 
  b.sec = 2*mu/(2*sigma^2)
  
  a.tert = ((3*mu)^2)/(3*sigma^2) 
  b.tert = 3*mu/(3*sigma^2)
  
  # coprimary transmission: both infected by the same primary case. |x2-x1| -> folded gamma difference dist
  
  f10 <- function(x, a.pri, b.pri)    (2-2*x)*dfgd(x, a = a.pri, b = b.pri)
  f1lower <- function(x,r,a.pri, b.pri)    	(x-r+1)*dfgd(x, a = a.pri, b = b.pri)
  f1upper <- function(x,r,a.pri, b.pri)    	(r+1-x)*dfgd(x, a = a.pri, b = b.pri)
  
  
  # primary - secondary. here it follows a gamma distribution 
  f20<-function(x,a.pri, b.pri)      	(2-2*x)*dgamma(x, shape = a.pri, rate = b.pri)
  f2lower<-function(x,r,a.pri, b.pri) 	(x-r+1)*dgamma(x, shape = a.pri, rate = b.pri)
  f2upper<-function(x,r,a.pri, b.pri) 	(r+1-x)*dgamma(x, shape = a.pri, rate = b.pri)
  
  f30<-function(x, a.sec, b.sec)       	(2-2*x)*dgamma(x, shape = a.sec, rate = b.sec)
  f3lower<-function(x,r,a.sec, b.sec) 	(x-r+1)*dgamma(x, shape = a.sec, rate = b.sec)
  f3upper<-function(x,r,a.sec, b.sec) 	(r+1-x)*dgamma(x, shape = a.sec, rate = b.sec)
  
  f40<-function(x, a.tert, b.tert)       	(2-2*x)*dgamma(x, shape = a.tert, rate = b.tert)
  f4lower<-function(x,r,a.tert, b.tert) 	(x-r+1)*dgamma(x, shape = a.tert, rate = b.tert)
  f4upper<-function(x,r,a.tert, b.tert) 	(r+1-x)*dgamma(x, shape = a.tert, rate = b.tert)
  
  #discretization
  p10 <- function(d, a.pri, b.pri)        	integrate(f = f10, lower = d, upper = (d+1), a.pri = a.pri, b.pri = b.pri)
  p1lower <- function(d, a.pri, b.pri)    	integrate(f = f1lower, lower = (d-1), upper = d, r = d, a.pri = a.pri, b.pri = b.pri)
  p1upper <- function(d, a.pri, b.pri)   	  integrate(f = f1upper, lower = d, upper = (d+1), r = d, a.pri = a.pri, b.pri = b.pri)
  
  p20 <- function(d, a.pri, b.pri)     	integrate(f = f20, lower = d, upper = (d+1), a.pri = a.pri, b.pri = b.pri)
  p2lower <- function(d, a.pri, b.pri) 	integrate(f = f2lower, lower = (d-1), upper = d, r = d, a.pri = a.pri, b.pri = b.pri)
  p2upper <- function(d, a.pri, b.pri) 	integrate(f = f2upper, lower = d, upper = (d+1), r = d, a.pri = a.pri, b.pri = b.pri)
  
  p30 <- function(d, a.sec, b.sec)     	integrate(f = f30, lower = d, upper = (d+1), a.sec = a.sec, b.sec = b.sec)
  p3lower <- function(d, a.sec, b.sec) 	integrate(f = f3lower, lower = (d-1), upper = d, r = d, a.sec = a.sec, b.sec = b.sec)
  p3upper <- function(d, a.sec, b.sec) 	integrate(f = f3upper, lower = d, upper = (d+1), r = d, a.sec = a.sec, b.sec = b.sec)
  
  p40 <- function(d, a.tert, b.tert)    	integrate(f = f40, lower = d, upper = (d+1), a.tert = a.tert, b.tert = b.tert)
  p4lower <- function(d, a.tert, b.tert) 	integrate(f=f4lower, lower=(d-1), upper=d, r = d, a.tert = a.tert, b.tert = b.tert)
  p4upper <- function(d, a.tert, b.tert) 	integrate(f=f4upper, lower=d, upper=(d+1), r = d, a.tert = a.tert, b.tert = b.tert)
  
  if(x==0){
    d1 <- p10(x, a.pri, b.pri)[[1]]
    d2 <- p20(x, a.pri, b.pri)[[1]]
    d3 <- p30(x, a.sec, b.sec)[[1]]
    d4 <- p40(x, a.tert, b.tert)[[1]]
  }
  else{
    d1<-p1lower(x, a.pri, b.pri)[[1]] + p1upper(x, a.pri, b.pri)[[1]]
    d2<-p2lower(x, a.pri, b.pri)[[1]] + p2upper(x, a.pri, b.pri)[[1]]
    d3<-p3lower(x, a.sec, b.sec)[[1]] + p3upper(x, a.sec, b.sec)[[1]]
    d4<-p4lower(x, a.tert, b.tert)[[1]] + p4upper(x, a.tert, b.tert)[[1]]
  }
  
  ind.output <- which.max(c(d1*w1, d2*w2, d3*w3, d4*w4))
  
  return(ind.output)
  
}

check_likelihood(20, GI_est_household[[100]]$est)
test_vector <- seq(3, 10, 1)
lapply(test_vector, function(x) check_likelihood(x, GI_est_household[[100]]$est))
### seems to work

### then for each sampling, the ICC intervals, check likelihood for their
### corresponding transmission generation
### 1 for co-prim, 2 for prim-sec, 3 for prim-tert, 4 for prim-quad

assign_cluster <- function(inferred_times, GIest){
  
  times.sort <- sort(inferred_times)
  n_case <- length(times.sort)
  ICCI <- times.sort[2:n_case] - times.sort[1]
  likelihood_checks <- unlist(lapply(
    ICCI, function(x) check_likelihood(x, GIest)))
  
  route1 <- length(which(likelihood_checks == 1))
  route2 <- length(which(likelihood_checks == 2))
  route3 <- length(which(likelihood_checks == 3))
  route4 <- length(which(likelihood_checks == 4))
  return(c(route1, route2, route3, route4))
}

assign_cluster(sample_times_list$household[[1]]$inferred.Tinfection[[1]],
               GI_est_household[[1]]$est)

View(sample_times_list)

length(sample_times_list$household[[1]]$inferred.Tinfection)

household_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 1782, ncol = 4)
  for(j in 1:1782){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$household[[i]]$inferred.Tinfection[[j]],
      GI_est_household[[i]]$est
    )
  }
  household_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}


length(sample_times_list$office[[1]]$inferred.Tinfection)

load("GI_est_office.rda")

office_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 88, ncol = 4)
  for(j in 1:88){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$office[[i]]$inferred.Tinfection[[j]],
      GI_est_office[[i]]$est
    )
  }
  office_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}

View(office_offspring_clust[[1]])
colMeans(office_offspring_clust[[1]])


load("GI_est_restaurant.rda")
length(sample_times_list$restaurant[[1]]$inferred.Tinfection)
# 175
restaurant_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 175, ncol = 4)
  for(j in 1:175){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$restaurant[[i]]$inferred.Tinfection[[j]],
      GI_est_restaurant[[i]]$est
    )
  }
  restaurant_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}

View(restaurant_offspring_clust[[1]])
colMeans(restaurant_offspring_clust[[2]])


load("GI_est_manual.rda")
length(sample_times_list$manual_labour[[1]]$inferred.Tinfection)
# 75
manual_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 75, ncol = 4)
  for(j in 1:75){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$manual_labour[[i]]$inferred.Tinfection[[j]],
      GI_est_manual[[i]]$est
    )
  }
  manual_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}

View(manual_offspring_clust[[1]])
colMeans(manual_offspring_clust[[1]])


load("GI_est_retail.rda")
length(sample_times_list$retail[[1]]$inferred.Tinfection)
# 17
retail_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 17, ncol = 4)
  for(j in 1:17){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$retail[[i]]$inferred.Tinfection[[j]],
      GI_est_retail[[i]]$est
    )
  }
  retail_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}


load("GI_est_noso.rda")
length(sample_times_list$noso[[1]]$inferred.Tinfection)
# 34
noso_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 34, ncol = 4)
  for(j in 1:34){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$noso[[i]]$inferred.Tinfection[[j]],
      GI_est_noso[[i]]$est
    )
  }
  noso_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}

View(noso_offspring_clust[[1]])


load("GI_est_social.rda")
length(sample_times_list$social[[1]]$inferred.Tinfection)
# 34
social_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 34, ncol = 4)
  for(j in 1:34){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$social[[i]]$inferred.Tinfection[[j]],
      GI_est_social[[i]]$est
    )
  }
  social_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}

View(social_offspring_clust[[1]])


load("GI_est_care.rda")
length(sample_times_list$carehome[[1]]$inferred.Tinfection)
# 23
carehome_offspring_clust <- vector("list", 100)

for(i in 1:100){
  offspring_matrix <- matrix(nrow = 23, ncol = 4)
  for(j in 1:23){
    offspring_matrix[j,] <- assign_cluster(
      sample_times_list$carehome[[i]]$inferred.Tinfection[[j]],
      GI_est_care[[i]]$est
    )
  }
  carehome_offspring_clust[[i]] <- offspring_matrix 
  print(i)
  rm(offspring_matrix)
}

offspring_infer_list <- list(
  household = household_offspring_clust,
  office = office_offspring_clust,
  restaurant = restaurant_offspring_clust,
  manual_labour = manual_offspring_clust,
  retail = retail_offspring_clust,
  nosocomial = noso_offspring_clust,
  social = social_offspring_clust,
  carehome = carehome_offspring_clust
)

save(offspring_infer_list, file = "offspring_infer_list.rda")

load("offspring_infer_list.rda")
####### now apply the model functions to estimate k and R in each setting
###### under a function that first deal with these clusters

source("tupper.R")

outputres <- function(modelres){
  data.frame(
    dispersion_k = modelres$a_n_b.est[1],
    meanclust = 1 + modelres$a_n_b.est[1] * modelres$a_n_b.est[2],
    Rc = modelres$a_n_b.est[1] * modelres$a_n_b.est[2]
  )
}


est_k_n_mu <- function(inputdata, n_single, q1, q2){
  offsprings <- as.vector(inputdata)
  offsprings <- offsprings[offsprings != 0]
  mod <- optim_fit_model(init = c(0.5, 0.5), 
                         inputdata = list(q1 = q1, 
                                          q2 = q2, 
                                          clusters = c(offsprings+1,
                                                       rep(1, n_single))))
  result <- outputres(mod)
  return(data.frame(
    k = as.numeric(result[1]),
    R = as.numeric(result[3])
  ))
}

est_k_n_mu(inputdata = household_offspring_clust[[20]], n_single = 1536,
           q1 = 1, q2 = 0.9)


household_offspring_clust <- offspring_infer_list$household
length(household_offspring_clust[[2]])


household_est <- vector("list", 100)
household_est_boot <- vector("list", 100)

for(i in 1:100){
  household_est[[i]] <- est_k_n_mu(household_offspring_clust[[i]],
                                   n_single = 1536,
                                   q1 = 1, q2 = 0.9)
  print(i)
}

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- household_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                n_single = 1536,
                                q1 = 1, q2 = 0.9), error = function(e) NA)
  }
  household_est_boot[[i]] <- est_boot
  print(i)
}

household_est_boot[[1]][[1]]

household_est_boot_k <- vector("list", 100)
for(i in 1:100){
  household_est_boot_k[[i]] <- unlist(lapply(household_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(household_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
#   18.01172 279.54098 

household_est_boot_R <- vector("list", 100)
for(i in 1:100){
  household_est_boot_R[[i]] <- unlist(lapply(household_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(household_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
#   0.7895574 0.8376054 

household_k <- lapply(household_est, function(x) x[1])
mean(unlist(household_k)) 
# 126.1795

household_R <- lapply(household_est, function(x) x[2])
mean(unlist(household_R)) 
# 0.8130676

household_offspring_clust <- offspring_infer_list$household

household_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- household_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 1536,
                                          q1 = 0.9, q2 = 1), error = function(e) NA)
  }
  household_est_boot_alt1[[i]] <- est_boot
  print(i)
}

household_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  household_est_boot_k_alt1[[i]] <- unlist(lapply(household_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(household_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
#   35.70647 390.16310 

household_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  household_est_boot_R_alt1[[i]] <- unlist(lapply(household_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(household_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.9401891 0.9967200 

household_est_alt1 <- vector("list", 100)

for(i in 1:100){
  household_est_alt1[[i]] <- est_k_n_mu(household_offspring_clust[[i]],
                                   n_single = 1536,
                                   q1 = 0.9, q2 = 1)
  print(i)
}

household_k_alt1 <- lapply(household_est_alt1, function(x) x[1])
mean(unlist(household_k_alt1))
# 185.3439

household_R_alt1 <- lapply(household_est_alt1, function(x) x[2])
mean(unlist(household_R_alt1))
# 0.9677

household_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  household_est_boot_R_alt1[[i]] <- unlist(lapply(household_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(household_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.9401891 0.9967200 

household_est_alt1 <- vector("list", 100)

for(i in 1:100){
  household_est_alt1[[i]] <- est_k_n_mu(household_offspring_clust[[i]],
                                        n_single = 1536,
                                        q1 = 0.9, q2 = 1)
  print(i)
}

household_k_alt1 <- lapply(household_est_alt1, function(x) x[1])
mean(unlist(household_k_alt1))
# 185.3439

household_R_alt1 <- lapply(household_est_alt1, function(x) x[2])
mean(unlist(household_R_alt1))
# 0.9677

household_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- household_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 1536,
                                          q1 = 0.9, q2 = 0.9), error = function(e) NA)
  }
  household_est_boot_alt2[[i]] <- est_boot
  print(i)
}

household_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  household_est_boot_k_alt2[[i]] <- unlist(lapply(household_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(household_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
#   19.93692 341.62148 

household_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  household_est_boot_R_alt2[[i]] <- unlist(lapply(household_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(household_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.8979456 0.9530981 

household_est_alt2 <- vector("list", 100)

for(i in 1:100){
  household_est_alt2[[i]] <- est_k_n_mu(household_offspring_clust[[i]],
                                        n_single = 1536,
                                        q1 = 0.9, q2 = 0.9)
  print(i)
}

household_k_alt2 <- lapply(household_est_alt2, function(x) x[1])
mean(unlist(household_k_alt2))
# 142.0361

household_R_alt2 <- lapply(household_est_alt2, function(x) x[2])
mean(unlist(household_R_alt2))
# 0.9249894

restaurant_offspring_clust <- offspring_infer_list$restaurant
restaurant_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- restaurant_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 107,
                                          q1 = 1, q2 = 0.8), error = function(e) NA)
  }
  restaurant_est_boot[[i]] <- est_boot
  print(i)
}

restaurant_est_boot_k <- vector("list", 100)
for(i in 1:100){
  restaurant_est_boot_k[[i]] <- unlist(lapply(restaurant_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(restaurant_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.9234479 3.3975088 

restaurant_est_boot_R <- vector("list", 100)
for(i in 1:100){
  restaurant_est_boot_R[[i]] <- unlist(lapply(restaurant_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(restaurant_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.9783081 1.3687812

restaurant_est <- vector("list", 100)

for(i in 1:100){
  restaurant_est[[i]] <- est_k_n_mu(restaurant_offspring_clust[[i]],
                                   n_single = 107,
                                   q1 = 1, q2 = 0.8)
  print(i)
}

restaurant_k <- lapply(restaurant_est, function(x) x[1])
mean(unlist(restaurant_k))
# 1.350834

restaurant_R <- lapply(restaurant_est, function(x) x[2])
mean(unlist(restaurant_R))
# 1.176404

restaurant_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- restaurant_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 107,
                                          q1 = 0.8, q2 = 1), error = function(e) NA)
  }
  restaurant_est_boot_alt1[[i]] <- est_boot
  print(i)
}

restaurant_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  restaurant_est_boot_k_alt1[[i]] <- unlist(lapply(restaurant_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(restaurant_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 1.054188 4.489080 

restaurant_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  restaurant_est_boot_R_alt1[[i]] <- unlist(lapply(restaurant_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(restaurant_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 1.388371 1.883866 

restaurant_est_alt1 <- vector("list", 100)

for(i in 1:100){
  restaurant_est_alt1[[i]] <- est_k_n_mu(restaurant_offspring_clust[[i]],
                                    n_single = 107,
                                    q1 = 0.8, q2 = 1)
  print(i)
}

restaurant_k_alt1 <- lapply(restaurant_est_alt1, function(x) x[1])
mean(unlist(restaurant_k_alt1))
# 1.582056

restaurant_R_alt1 <- lapply(restaurant_est_alt1, function(x) x[2])
mean(unlist(restaurant_R_alt1))
# 1.639793

restaurant_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- restaurant_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 107,
                                          q1 = 0.8, q2 = 0.8), error = function(e) NA)
  }
  restaurant_est_boot_alt2[[i]] <- est_boot
  print(i)
}

restaurant_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  restaurant_est_boot_k_alt2[[i]] <- unlist(lapply(restaurant_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(restaurant_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.9095643 3.4851476 

restaurant_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  restaurant_est_boot_R_alt2[[i]] <- unlist(lapply(restaurant_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(restaurant_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 1.26360 1.70697 

restaurant_est_alt2 <- vector("list", 100)

for(i in 1:100){
  restaurant_est_alt2[[i]] <- est_k_n_mu(restaurant_offspring_clust[[i]],
                                         n_single = 107,
                                         q1 = 0.8, q2 = 0.8)
  print(i)
}

restaurant_k_alt2 <- lapply(restaurant_est_alt2, function(x) x[1])
mean(unlist(restaurant_k_alt2))
# 1.346139

restaurant_R_alt2 <- lapply(restaurant_est_alt2, function(x) x[2])
mean(unlist(restaurant_R_alt2))
# 1.4868


carehome_offspring_clust <- offspring_infer_list$carehome
carehome_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- carehome_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 36,
                                          q1 = 1, q2 = 0.7), error = function(e) NA)
  }
  carehome_est_boot[[i]] <- est_boot
  print(i)
}

carehome_est_boot_k <- vector("list", 100)
for(i in 1:100){
  carehome_est_boot_k[[i]] <- unlist(lapply(carehome_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(carehome_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.1898167 0.4094464 

carehome_est_boot_R <- vector("list", 100)
for(i in 1:100){
  carehome_est_boot_R[[i]] <- unlist(lapply(carehome_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(carehome_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 1.633801 4.752922 

carehome_est <- vector("list", 100)

for(i in 1:100){
  carehome_est[[i]] <- est_k_n_mu(carehome_offspring_clust[[i]],
                                       n_single = 36,
                                       q1 = 1, q2 = 0.7)
  print(i)
}

carehome_k <- lapply(carehome_est, function(x) x[1])
mean(unlist(carehome_k))
# 0.2716658

carehome_R <- lapply(carehome_est, function(x) x[2])
mean(unlist(carehome_R))
# 2.637834

carehome_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- carehome_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 36,
                                          q1 = 0.9, q2 = 1), error = function(e) NA)
  }
  carehome_est_boot_alt1[[i]] <- est_boot
  print(i)
}

carehome_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  carehome_est_boot_k_alt1[[i]] <- unlist(lapply(carehome_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(carehome_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.2336222 0.5113960 

carehome_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  carehome_est_boot_R_alt1[[i]] <- unlist(lapply(carehome_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(carehome_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 2.134699 5.932875 

carehome_est_alt1 <- vector("list", 100)

for(i in 1:100){
  carehome_est_alt1[[i]] <- est_k_n_mu(carehome_offspring_clust[[i]],
                                    n_single = 36,
                                    q1 = 0.9, q2 = 1)
  print(i)
}

carehome_k_alt1 <- lapply(carehome_est_alt1, function(x) x[1])
mean(unlist(carehome_k_alt1))
# 0.3346867

carehome_R_alt1 <- lapply(carehome_est_alt1, function(x) x[2])
mean(unlist(carehome_R_alt1))
# 3.40

carehome_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- carehome_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 36,
                                          q1 = 0.9, q2 = 0.7), error = function(e) NA)
  }
  carehome_est_boot_alt2[[i]] <- est_boot
  print(i)
}

carehome_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  carehome_est_boot_k_alt2[[i]] <- unlist(lapply(carehome_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(carehome_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.1820572 0.3992437 

carehome_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  carehome_est_boot_R_alt2[[i]] <- unlist(lapply(carehome_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(carehome_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 1.767788 5.025942 

carehome_est_alt2 <- vector("list", 100)

for(i in 1:100){
  carehome_est_alt2[[i]] <- est_k_n_mu(carehome_offspring_clust[[i]],
                                       n_single = 36,
                                       q1 = 0.9, q2 = 0.7)
  print(i)
}

carehome_k_alt2 <- lapply(carehome_est_alt2, function(x) x[1])
mean(unlist(carehome_k_alt2))
# 0.2623302

carehome_R_alt2 <- lapply(carehome_est_alt2, function(x) x[2])
mean(unlist(carehome_R_alt2))
# 2.822933

noso_offspring_clust <- offspring_infer_list$nosocomial

noso_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- noso_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 46,
                                          q1 = 1, q2 = 0.9), error = function(e) NA)
  }
  noso_est_boot[[i]] <- est_boot
  print(i)
}

noso_est_boot_k <- vector("list", 100)
for(i in 1:100){
  noso_est_boot_k[[i]] <- unlist(lapply(noso_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(noso_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.3973615 3.4079625 

noso_est_boot_R <- vector("list", 100)
for(i in 1:100){
  noso_est_boot_R[[i]] <- unlist(lapply(noso_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(noso_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.725913 1.469557


nosocomial_est <- vector("list", 100)

for(i in 1:100){
  nosocomial_est[[i]] <- est_k_n_mu(noso_offspring_clust[[i]],
                                  n_single = 46,
                                  q1 = 1, q2 = 0.9)
  print(i)
}

nosocomial_k <- lapply(nosocomial_est, function(x) x[1])
mean(unlist(nosocomial_k))
# 0.7605366

nosocomial_R <- lapply(nosocomial_est, function(x) x[2])
mean(unlist(nosocomial_R))
# 1.03

noso_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- noso_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 46,
                                          q1 = 0.9, q2 = 1), error = function(e) NA)
  }
  noso_est_boot_alt1[[i]] <- est_boot
  print(i)
}

noso_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  noso_est_boot_k_alt1[[i]] <- unlist(lapply(noso_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(noso_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.4230654 4.4450056 

noso_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  noso_est_boot_R_alt1[[i]] <- unlist(lapply(noso_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(noso_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.8629946 1.6837906 

nosocomial_est_alt1 <- vector("list", 100)

for(i in 1:100){
  nosocomial_est_alt1[[i]] <- est_k_n_mu(noso_offspring_clust[[i]],
                                    n_single = 46,
                                    q1 = 0.9, q2 = 1)
  print(i)
}

nosocomial_k_alt1 <- lapply(nosocomial_est_alt1, function(x) x[1])
mean(unlist(nosocomial_k_alt1))
# 0.8270169

nosocomial_R_alt1 <- lapply(nosocomial_est_alt1, function(x) x[2])
mean(unlist(nosocomial_R_alt1))
# 1.206082

noso_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- noso_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 46,
                                          q1 = 0.9, q2 = 0.9), error = function(e) NA)
  }
  noso_est_boot_alt2[[i]] <- est_boot
  print(i)
}

noso_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  noso_est_boot_k_alt2[[i]] <- unlist(lapply(noso_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(noso_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.3882808 3.5300248 

noso_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  noso_est_boot_R_alt2[[i]] <- unlist(lapply(noso_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(noso_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.8157403 1.5955860 

nosocomial_est_alt2 <- vector("list", 100)

for(i in 1:100){
  nosocomial_est_alt2[[i]] <- est_k_n_mu(noso_offspring_clust[[i]],
                                         n_single = 46,
                                         q1 = 0.9, q2 = 0.9)
  print(i)
}

nosocomial_k_alt2 <- lapply(nosocomial_est_alt2, function(x) x[1])
mean(unlist(nosocomial_k_alt2))
# 0.7542489

nosocomial_R_alt2 <- lapply(nosocomial_est_alt2, function(x) x[2])
mean(unlist(nosocomial_R_alt2))
# 1.14

######### TBD 20250311

manual_offspring_clust <- offspring_infer_list$manual_labour

manual_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- manual_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 178,
                                          q1 = 1, q2 = 0.8), error = function(e) NA)
  }
  manual_est_boot[[i]] <- est_boot
  print(i)
}

manual_est_boot_k <- vector("list", 100)
for(i in 1:100){
  manual_est_boot_k[[i]] <- unlist(lapply(manual_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(manual_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.1775539 0.4945517 

manual_est_boot_R <- vector("list", 100)
for(i in 1:100){
  manual_est_boot_R[[i]] <- unlist(lapply(manual_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(manual_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.6272783 1.3078441 

manual_est <- vector("list", 100)

for(i in 1:100){
  manual_est[[i]] <- est_k_n_mu(manual_offspring_clust[[i]],
                                    n_single = 178,
                                    q1 = 1, q2 = 0.8)
  print(i)
}

manual_k <- lapply(manual_est, function(x) x[1])
mean(unlist(manual_k))
# 0.2544849

manual_R <- lapply(manual_est, function(x) x[2])
mean(unlist(manual_R))
# 0.92

manual_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- manual_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 178,
                                          q1 = 0.8, q2 = 1), error = function(e) NA)
  }
  manual_est_boot_alt1[[i]] <- est_boot
  print(i)
}

manual_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  manual_est_boot_k_alt1[[i]] <- unlist(lapply(manual_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(manual_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.1940823 0.5779232 

manual_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  manual_est_boot_R_alt1[[i]] <- unlist(lapply(manual_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(manual_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.8748885 1.7293688 

manual_est_alt1 <- vector("list", 100)

for(i in 1:100){
  manual_est_alt1[[i]] <- est_k_n_mu(manual_offspring_clust[[i]],
                                n_single = 178,
                                q1 = 0.8, q2 = 1)
  print(i)
}

manual_k_alt1 <- lapply(manual_est_alt1, function(x) x[1])
mean(unlist(manual_k_alt1))
# 0.2834124

manual_R_alt1 <- lapply(manual_est_alt1, function(x) x[2])
mean(unlist(manual_R_alt1))
# 1.240506

manual_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- manual_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 178,
                                          q1 = 0.8, q2 = 0.8), error = function(e) NA)
  }
  manual_est_boot_alt2[[i]] <- est_boot
  print(i)
}

manual_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  manual_est_boot_k_alt2[[i]] <- unlist(lapply(manual_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(manual_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.1639898 0.4798888 

manual_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  manual_est_boot_R_alt2[[i]] <- unlist(lapply(manual_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(manual_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.7534266 1.4828723 

manual_est_alt2 <- vector("list", 100)

for(i in 1:100){
  manual_est_alt2[[i]] <- est_k_n_mu(manual_offspring_clust[[i]],
                                     n_single = 178,
                                     q1 = 0.8, q2 = 0.8)
  print(i)
}

manual_k_alt2 <- lapply(manual_est_alt2, function(x) x[1])
mean(unlist(manual_k_alt2))
# 0.239485

manual_R_alt2 <- lapply(manual_est_alt2, function(x) x[2])
mean(unlist(manual_R_alt2))
# 1.063943


office_offspring_clust <- offspring_infer_list$office

office_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- office_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 277,
                                          q1 = 1, q2 = 0.9), error = function(e) NA)
  }
  office_est_boot[[i]] <- est_boot
  print(i)
}

office_est_boot_k <- vector("list", 100)
for(i in 1:100){
  office_est_boot_k[[i]] <- unlist(lapply(office_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(office_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.8405767 123.7576611 

office_est_boot_R <- vector("list", 100)
for(i in 1:100){
  office_est_boot_R[[i]] <- unlist(lapply(office_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(office_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.2954946 0.3972215

office_est <- vector("list", 100)

for(i in 1:100){
  office_est[[i]] <- est_k_n_mu(office_offspring_clust[[i]],
                                n_single = 277,
                                q1 = 1, q2 = 0.9)
  print(i)
}

office_k <- lapply(office_est, function(x) x[1])
mean(unlist(office_k))
# 9.117365

office_R <- lapply(office_est, function(x) x[2])
mean(unlist(office_R))
# 0.34044

office_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- office_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 277,
                                          q1 = 0.9, q2 = 1), error = function(e) NA)
  }
  office_est_boot_alt1[[i]] <- est_boot
  print(i)
}

office_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  office_est_boot_k_alt1[[i]] <- unlist(lapply(office_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(office_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
#0.9707206 149.7881953 

office_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  office_est_boot_R_alt1[[i]] <- unlist(lapply(office_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(office_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.3570117 0.4750495 

office_est_alt1 <- vector("list", 100)

for(i in 1:100){
  office_est_alt1[[i]] <- est_k_n_mu(office_offspring_clust[[i]],
                                n_single = 277,
                                q1 = 0.9, q2 = 1)
  print(i)
}

office_k_alt1 <- lapply(office_est_alt1, function(x) x[1])
mean(unlist(office_k_alt1))
# 13.52804

office_R_alt1 <- lapply(office_est_alt1, function(x) x[2])
mean(unlist(office_R_alt1))
# 0.41028

office_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- office_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 277,
                                          q1 = 0.9, q2 = 0.9), error = function(e) NA)
  }
  office_est_boot_alt2[[i]] <- est_boot
  print(i)
}

office_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  office_est_boot_k_alt2[[i]] <- unlist(lapply(office_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(office_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
#0.8397473 134.0069511

office_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  office_est_boot_R_alt2[[i]] <- unlist(lapply(office_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(office_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.3309975 0.4427960 

office_est_alt2 <- vector("list", 100)

for(i in 1:100){
  office_est_alt2[[i]] <- est_k_n_mu(office_offspring_clust[[i]],
                                     n_single = 277,
                                     q1 = 0.9, q2 = 0.9)
  print(i)
}

office_k_alt2 <- lapply(office_est_alt2, function(x) x[1])
mean(unlist(office_k_alt2))
# 9.648378

office_R_alt2 <- lapply(office_est_alt2, function(x) x[2])
mean(unlist(office_R_alt2))
# 0.3817333


social_offspring_clust <- offspring_infer_list$social

social_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- social_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 27,
                                          q1 = 1, q2 = 0.5), error = function(e) NA)
  }
  social_est_boot[[i]] <- est_boot
  print(i)
}

social_est_boot_k <- vector("list", 100)
for(i in 1:100){
  social_est_boot_k[[i]] <- unlist(lapply(social_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(social_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.1249585 0.4661757 

social_est_boot_R <- vector("list", 100)
for(i in 1:100){
  social_est_boot_R[[i]] <- unlist(lapply(social_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(social_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 1.373492 11.666623

social_est <- vector("list", 100)

for(i in 1:100){
  social_est[[i]] <- est_k_n_mu(social_offspring_clust[[i]],
                                n_single = 27,
                                q1 = 1, q2 = 0.5)
  print(i)
}

social_k <- lapply(social_est, function(x) x[1])
mean(unlist(social_k))
# 0.1670893

social_R <- lapply(social_est, function(x) x[2])
mean(unlist(social_R))
# 5.32096

social_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- social_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 27,
                                          q1 = 0.7, q2 = 1), error = function(e) NA)
  }
  social_est_boot_alt1[[i]] <- est_boot
  print(i)
}

social_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  social_est_boot_k_alt1[[i]] <- unlist(lapply(social_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(social_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 0.1590934 0.6203980 

social_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  social_est_boot_R_alt1[[i]] <- unlist(lapply(social_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(social_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.625136 20.724352 

social_est_alt1 <- vector("list", 100)

for(i in 1:100){
  social_est_alt1[[i]] <- est_k_n_mu(social_offspring_clust[[i]],
                                n_single = 27,
                                q1 = 0.7, q2 = 1)
  print(i)
}

social_k_alt1 <- lapply(social_est_alt1, function(x) x[1])
mean(unlist(social_k_alt1))
# 0.2138631

social_R_alt1 <- lapply(social_est_alt1, function(x) x[2])
mean(unlist(social_R_alt1))
# 9.606909


social_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- social_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 27,
                                          q1 = 0.7, q2 = 0.5), error = function(e) NA)
  }
  social_est_boot_alt2[[i]] <- est_boot
  print(i)
}

social_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  social_est_boot_k_alt2[[i]] <- unlist(lapply(social_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(social_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 

social_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  social_est_boot_R_alt2[[i]] <- unlist(lapply(social_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(social_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 

social_est_alt2 <- vector("list", 100)

for(i in 1:100){
  social_est_alt2[[i]] <- est_k_n_mu(social_offspring_clust[[i]],
                                     n_single = 27,
                                     q1 = 0.7, q2 = 0.5)
  print(i)
}

social_k_alt2 <- lapply(social_est_alt2, function(x) x[1])
mean(unlist(social_k_alt2))
# 2.5%      50%    97.5% 
# 0.2019465 0.2145959 0.2245658

social_R_alt2 <- lapply(social_est_alt2, function(x) x[2])
mean(unlist(social_R_alt2))
# 2.5%      50%    97.5% 
# 8.907804  9.608734 10.276642 


retail_offspring_clust <- offspring_infer_list$retail
  
retail_est_boot <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- retail_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 91,
                                          q1 = 1, q2 = 0.7), error = function(e) NA)
  }
  retail_est_boot[[i]] <- est_boot
  print(i)
}

retail_est_boot_k <- vector("list", 100)
for(i in 1:100){
  retail_est_boot_k[[i]] <- unlist(lapply(retail_est_boot[[i]], function(x) x[1]))
}
quantile(unlist(retail_est_boot_k), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.05257211 0.34644481 

retail_est_boot_R <- vector("list", 100)
for(i in 1:100){
  retail_est_boot_R[[i]] <- unlist(lapply(retail_est_boot[[i]], function(x) x[2]))
}
quantile(unlist(retail_est_boot_R), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.2434879 1.8199286 


retail_est <- vector("list", 100)

for(i in 1:100){
  retail_est[[i]] <- est_k_n_mu(retail_offspring_clust[[i]],
                                n_single = 91,
                                q1 = 1, q2 = 0.7)
  print(i)
}

retail_k <- lapply(retail_est, function(x) x[1])
mean(unlist(retail_k))
# 0.1015475

retail_R <- lapply(retail_est, function(x) x[2])
mean(unlist(retail_R))
# 0.5497724

retail_est_boot_alt1 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- retail_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 91,
                                          q1 = 0.7, q2 = 1), error = function(e) NA)
  }
  retail_est_boot_alt1[[i]] <- est_boot
  print(i)
}

retail_est_boot_k_alt1 <- vector("list", 100)
for(i in 1:100){
  retail_est_boot_k_alt1[[i]] <- unlist(lapply(retail_est_boot_alt1[[i]], function(x) x[1]))
}
quantile(unlist(retail_est_boot_k_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.05798939 0.43292009 

retail_est_boot_R_alt1 <- vector("list", 100)
for(i in 1:100){
  retail_est_boot_R_alt1[[i]] <- unlist(lapply(retail_est_boot_alt1[[i]], function(x) x[2]))
}
quantile(unlist(retail_est_boot_R_alt1), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.4195563 2.5200156 


retail_est_alt1 <- vector("list", 100)

for(i in 1:100){
  retail_est_alt1[[i]] <- est_k_n_mu(retail_offspring_clust[[i]],
                                n_single = 91,
                                q1 = 0.7, q2 = 1)
  print(i)
}

retail_k_alt1 <- lapply(retail_est_alt1, function(x) x[1])
mean(unlist(retail_k_alt1))
# 0.1154391

retail_R_alt1 <- lapply(retail_est_alt1, function(x) x[2])
mean(unlist(retail_R_alt1))
# 0.8643133


retail_est_boot_alt2 <- vector("list", 100)

for(i in 1:100){
  est_boot <- vector("list", 100)
  offsprings <- retail_offspring_clust[[i]]
  for(j in 1:100){
    set.seed(j)
    offsprings_boot <- sample(offsprings, length(offsprings), replace = T) 
    est_boot[[j]] <- tryCatch (est_k_n_mu(offsprings_boot,
                                          n_single = 91,
                                          q1 = 0.7, q2 = 0.7), error = function(e) NA)
  }
  retail_est_boot_alt2[[i]] <- est_boot
  print(i)
}

retail_est_boot_k_alt2 <- vector("list", 100)
for(i in 1:100){
  retail_est_boot_k_alt2[[i]] <- unlist(lapply(retail_est_boot_alt2[[i]], function(x) x[1]))
}
quantile(unlist(retail_est_boot_k_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.04362133 0.33081667 

retail_est_boot_R_alt2 <- vector("list", 100)
for(i in 1:100){
  retail_est_boot_R_alt2[[i]] <- unlist(lapply(retail_est_boot_alt2[[i]], function(x) x[2]))
}
quantile(unlist(retail_est_boot_R_alt2), c(0.025, 0.975), na.rm = T)
# 2.5%     97.5% 
# 0.4195563 2.5200156 

retail_est_alt2 <- vector("list", 100)

for(i in 1:100){
  retail_est_alt2[[i]] <- est_k_n_mu(retail_offspring_clust[[i]],
                                     n_single = 91,
                                     q1 = 0.7, q2 = 0.7)
  print(i)
}

retail_k_alt2 <- lapply(retail_est_alt2, function(x) x[1])
mean(unlist(retail_k_alt2))
# 0.08808471
retail_R_alt2 <- lapply(retail_est_alt2, function(x) x[2])
mean(unlist(retail_R_alt2))
# 0.643061










k_mu_ests <- list(
  household = household_est,
  restaurant = restaurant_est,
  carehome = carehome_est,
  nosocomial = nosocomial_est,
  manual_labour = manual_est,
  office = office_est,
  social = social_est,
  retail = retail_est
)

save(k_mu_ests, file = "k_n_mu_est_ind.rda")




### estimate prop_20_80

propresponsible <- function(R0, k, prop = 0.8) {
  qm1 <- qnbinom(1 - prop, k + 1, mu = R0 * (k + 1) / k)
  remq <- 1 - prop - pnbinom(qm1 - 1, k + 1, mu = R0 * (k + 1) / k)
  remx <- remq / dnbinom(qm1, k + 1, mu = R0 * (k + 1) / k)
  q <- qm1 + 1
  1 - pnbinom(q - 1, k, mu = R0) - dnbinom(q, k, mu = R0) * remx
}



frac_est_CI.2080 <- function(Kestinput, Restinput){
  frac_all = numeric(100)
  for(i in 1:100){
    frac_all[i] <- propresponsible(R0 = as.numeric(Restinput[[i]]),
                                   k = as.numeric(Kestinput[[i]]), prop = 0.8)
  }
  
  return(quantile(frac_all, c(0.025, 0.5, 0.975)))
}



frac_est_CI.2080(Kestinput = household_k,
                 Restinput = household_R)


# 2.5%       50%     97.5% 
# 0.3881695 0.3922751 0.3947944 


frac_est_CI.2080(Kestinput = restaurant_k,
                 Restinput = restaurant_R)

# 2.5%       50%     97.5% 
# 0.3211268 0.3352870 0.3462434 


frac_est_CI.2080(Kestinput = carehome_k,
                 Restinput = carehome_R)

# 2.5%       50%     97.5% 
# 0.1771906 0.1961145 0.2114462 

frac_est_CI.2080(Kestinput = nosocomial_k,
                 Restinput = nosocomial_R)
# 2.5%       50%     97.5% 
# 0.2445946 0.2734811 0.3028370 

frac_est_CI.2080(Kestinput = manual_k,
                 Restinput = manual_R)
# 2.5%       50%     97.5% 
# 0.1539109 0.1628534 0.1714555 

frac_est_CI.2080(Kestinput = office_k,
                 Restinput = office_R)
# 2.5%       50%     97.5% 
# 0.1950328 0.2077258 0.2173406 

frac_est_CI.2080(Kestinput = social_k,
                 Restinput = social_R)
# 2.5%       50%     97.5% 
# 0.1379062 0.1450707 0.1504519 

frac_est_CI.2080(Kestinput = retail_k,
                 Restinput = retail_R)
# 2.5%       50%     97.5% 
# 0.06673620 0.08201863 0.09259514 









