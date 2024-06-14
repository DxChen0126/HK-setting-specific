rm(list = ls())

# simulation code, derived from Park's code

##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
  if(length(x)==1) {
    rep(x, size)
  } else{
    sample(x, size, replace=TRUE)
  }
}

sir.full <- function(size,
                     R, # input R
                     k, # input dispersion k
                     para,
                     I0, # inial infectors
                     seed = NULL,
                     imax,
                     keep.intrinsic=FALSE,
                     GItrunc){
  # para is a vector contains of 
  # c(incumu, incusd, GImu, GIsd)
  
  # gamma shape and rate parameters of IP
  incushape = (para[1]/para[2])^2
  incurate = para[1]/para[2]^2
  
  # gamma shape and rate parameters of GI
  GIshape = (para[3]/para[4])^2
  GIrate = para[3]/para[4]^2
  
  if (!is.null(seed)) set.seed(seed)
  
  # initial population size
  V <- 1:size
  
  # specify the initial infectors
  if (missing(I0)) {
    if (missing(I0)) stop("specify the initial conditions")
    
  }
  initial_infected <- 1:I0
  
  # specify number of maximum infected
  if (missing(imax)) imax <- size
  
  queue_v <- queue_t <- queue_infector <- rep(NA, I0)
  
  queue_v[1:I0] <- initial_infected 
  queue_t[1:I0] <- 0
  
  t_infected <- t_symptomatic <- rep(NA, size)
  t_infected[initial_infected] <- 0
  
  t_gillespie <- NULL
  c_infected <- 0
  
  if (keep.intrinsic) {
    intrinsic_generation <- vector('list', length(V))
  } else {
    intrinsic_generation <- NULL
  }
  
  done <- rep(FALSE, size)
  infected_by <- rep(NA, size)
  
  stop <- FALSE
  
  while (!stop) {
    j.index <- which.min(queue_t)
    j <- queue_v[j.index]
    
    infected_by[j] <- queue_infector[j.index]
    t_infected[j] <- queue_t[j.index]
    
    t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
    
    # sample incubation period
    incubation <- rgamma(1, shape = incushape, rate = incurate)
    t_symptomatic[j] <- t+incubation
    
    c_infected <- c_infected +1
    
    # sample number of contacts
    # ncontact <- rpois(1, R0)
    # now use the negative binomial setting
    ncontact <- rnbinom(1, size = k, mu = R)
    
    # truncation of the generation interval also further impact ncontact
    
    #  sample generation interval
    generation <- rgamma(ncontact, shape = GIshape, rate = GIrate)
    
    # exclude cases beyond truncation threshold
    ind.exc <- which(generation > GItrunc)
    
    
    if(length(ind.exc) > 0){
    ncontact <- ncontact - length(ind.exc)
    generation <- generation[-ind.exc]
    }
    
    n <- V[V != j]
    
    if (ncontact > 0) {
      queue_v <- c(queue_v, sample2(n, ncontact))
      queue_infector <- c(queue_infector, rep(j, ncontact))
      queue_t <- c(queue_t, t + generation)
    }
    
    if (keep.intrinsic) intrinsic_generation[[j]] <- generation
    
    # if (ncontact > 0) {
    #   queue_t <- c(queue_t, t + generation)
    # }
    
    done[j] <- TRUE
    
    filter2 <- !done[queue_v]
    queue_v <- queue_v[filter2]
    queue_infector <- queue_infector[filter2]
    queue_t <- queue_t[filter2]
    
    stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected == imax)
  }
  
  return(
    list(
      data=data.frame(
        time=t_gillespie[(I0):c_infected],
        infected=(I0):c_infected
      ),
      intrinsic_generation=intrinsic_generation,
      t_infected=t_infected,
      t_symptomatic=t_symptomatic,
      infected_by=infected_by,
      maxsize = max(table(infected_by))
    )
  )
}

# parameter options

R_opt <- c(0.5, 1, 2)
k_opt <- c(0.1, 0.5, 1)
GItrunc_opt1 <- c(2, 4, 6)
GItrunc_opt2 <- c(8, 10, 12)


paracomb1 <- data.frame(
  R = rep(R_opt, each = 9),
  k = rep(rep(k_opt, each = 3), 3),
  GItrunc = rep(GItrunc_opt1, 9)
)

paracomb2 <- data.frame(
  R = rep(R_opt, each = 9),
  k = rep(rep(k_opt, each = 3), 3),
  GItrunc = rep(GItrunc_opt2, 9)
)


simreslist1 <- vector("list", 27)

progressbar <- txtProgressBar(min = 0, max = 27, style = 3)

for(j in 1:27){
  maxsize <- numeric(1000)
  for(i in 1:1000){
    simres <- sir.full(size = 1000,
                       R = paracomb1$R[j], # input R
                       k = paracomb1$k[j], # input dispersion k
                       para = c(6.5, 2.6, 7, 4),
                       I0 = 10, # initial infectors
                       seed = i,
                       keep.intrinsic=FALSE,
                       GItrunc = paracomb1$GItrunc[j])
    maxsize[i] <- simres$maxsize
  }
  ind.check <- which(is.infinite(maxsize))
  if(length(ind.check) > 0){maxsize[ind.check] <- 0}
  simreslist1[[j]] <- maxsize
  setTxtProgressBar(progressbar, j)
}

save(file = "sim_max_size_truncGI_veryshort.rda", simreslist1)

load("sim_max_size_truncGI_veryshort.rda")

simreslist2 <- vector("list", 27)

progressbar <- txtProgressBar(min = 0, max = 27, style = 3)

for(j in 1:27){
  maxsize <- numeric(1000)
  for(i in 1:1000){
    simres <- sir.full(size = 1000,
                       R = paracomb2$R[j], # input R
                       k = paracomb2$k[j], # input dispersion k
                       para = c(6.5, 2.6, 7, 4),
                       I0 = 10, # initial infectors
                       seed = i,
                       keep.intrinsic=FALSE,
                       GItrunc = paracomb2$GItrunc[j])
    maxsize[i] <- simres$maxsize
  }
  ind.check <- which(is.infinite(maxsize))
  if(length(ind.check) > 0){maxsize[ind.check] <- 0}
  simreslist2[[j]] <- maxsize
  setTxtProgressBar(progressbar, j)
}

save(file = "sim_max_size_truncGI_relativelong.rda", simreslist2)

load("sim_max_size_truncGI_relativelong.rda")

df1 <- data.frame(
  R = rep(rep(0.5, 9000), 2),
  k = rep(rep(k_opt, each = 3000), 2),
  GI_trunc = c(rep(rep(c(2, 4, 6), each = 1000), 3), rep(rep(c(8, 10, 12), each = 1000), 3)),
  maxsize = c(unlist(simreslist1[1:9]), unlist(simreslist2[1:9]))
)

df1$k <- factor(df1$k, levels = c(0.1, 0.5, 1))
df1$GI_trunc <- factor(df1$GI_trunc, levels = c(2, 4, 6, 8, 10, 12))



p1 <- ggplot(data = df1) +
  geom_violin(aes(fill = GI_trunc, y = maxsize, x = k), position = position_dodge(width = 0.75), scale="width") + 
  scale_fill_brewer(palette = "YlGnBu", name = "Truncated GI length") + 
  ylab("Maximum offspring size obtained") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) + 
  ggtitle("R = 0.5") 
p1

df2 <- data.frame(
  R = rep(rep(1, 9000), 2),
  k = rep(rep(k_opt, each = 3000), 2),
  GI_trunc = c(rep(rep(c(2, 4, 6), each = 1000), 3), rep(rep(c(8, 10, 12), each = 1000), 3)),
  maxsize = c(unlist(simreslist1[10:18]), unlist(simreslist2[10:18]))
)

df2$k <- factor(df2$k, levels = c(0.1, 0.5, 1))
df2$GI_trunc <- factor(df2$GI_trunc, levels = c(2, 4, 6, 8, 10, 12))



p2 <- ggplot(data = df2) +
  geom_violin(aes(fill = GI_trunc, y = maxsize, x = k), position = position_dodge(width = 0.75), scale="width") + 
  scale_fill_brewer(palette = "YlGnBu", name = "Truncated GI length") + 
  ylab("Maximum offspring size obtained") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) + 
  ggtitle("R = 1") 
p2

df3 <- data.frame(
  R = rep(rep(2, 9000), 2),
  k = rep(rep(k_opt, each = 3000), 2),
  GI_trunc = c(rep(rep(c(2, 4, 6), each = 1000), 3), rep(rep(c(8, 10, 12), each = 1000), 3)),
  maxsize = c(unlist(simreslist1[19:27]), unlist(simreslist2[19:27]))
)

df3$k <- factor(df3$k, levels = c(0.1, 0.5, 1))
df3$GI_trunc <- factor(df3$GI_trunc, levels = c(2, 4, 6, 8, 10, 12))



p3 <- ggplot(data = df3) +
  geom_violin(aes(fill = GI_trunc, y = maxsize, x = k), position = position_dodge(width = 0.75), scale="width") + 
  scale_fill_brewer(palette = "YlGnBu", name = "Truncated GI length") + 
  ylab("Maximum offspring size obtained") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) + 
  ggtitle("R = 2") 
p3

library(ggpubr)

psize_sim <- ggarrange(p1, p2, p3, nrow = 3, ncol = 1, common.legend = T)
ggsave("./event cluster analysis/slides/pfig_v4_slides.pdf", width = 12.5, height = 9, units = "in")
ggsave("maxsize_sim_20231115.pdf", width = 12.5, height = 9, units = "in")





