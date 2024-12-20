library(ggplot2)
library(ggpubr)


# density of gamma difference distribution
dgdd <- function(x,a,b){
  # here a and b are shape and rate
  if(length(a)!=2 || length(b)!=2) stop("Provide a=(a1,a2) and b=(b1,b2).")
  if(a[1]<=0 || a[2]<=0) stop("NaN produced for a<=0.")
  if(b[1]<=0 || b[2]<=0) stop("Provide b>0.")
  gt <- function(ti,a,b){
    cons <- b[1]^a[1] * b[2]^a[2] / (gamma(a[1])*gamma(a[2]))
    if(ti>=0){
      integrand <- function(x,s) x^(a[1]-1) * (x-s)^(a[2]-1) * exp(-sum(b)*x)
      int <- tryCatch(integrate(integrand, lower=ti, upper=Inf, s=ti)$value, error=function(e) return(0))
      #when the integral is divergent, set the density=0
      out <- cons*exp(b[2]*ti)*int
    } 
    else{
      integrand <- function(x,s) x^(a[2]-1) * (x+s)^(a[1]-1) * exp(-sum(b)*x)
      int <- tryCatch(integrate(integrand, lower=-ti, upper=Inf, s=ti)$value, error=function(e) return(0))
      #when the integral is divergent, set the density=0
      out <- cons*exp(-b[1]*ti)*int
    } 
    if(is.nan(out)) out <- 0 #happens because we may multiply 0 and Inf
    return(out)
  }
  ft <- tryCatch(sapply(1:length(x), function(i) gt(ti=x[i], a, b)),
                 error=function(e){
                   cat("Problem when evaluation x at:", x[i], "\n")
                   cat("print(par):", a, b, "\n")
                   msg <- conditionMessage(e)
                   return(msg)
                 })
  if(class(ft)!="numeric") stop("")
  return(ft)
}

# density of folded gamma difference
dfgd <- function(x,a,b){
  if(a<=0) stop("NaN produced for a<=0.")
  if(b<=0) stop("Provide b>0.")
  ft <- numeric(length(x))
  # note here c(a, a) c(b, b) means these two r.v.s are i.i.d from a same Gamma dist
  ft[which(x>=0)] <- dgdd(x[x>=0], c(a,a), c(b,b))
  return(2*ft)
}


GTest <- read_excel("update GI setting estimates.xlsx")


estimate <- GTest[4,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

load("sample_times_list.rda")
sample_times_list$household

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

lapply(sample_times_list$household, function(x) arrange_ICCI(x$inferred.Tinfection))

householdICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$household, function(x) arrange_ICCI(x$inferred.Tinfection)))
)


phousehold <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                       y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = householdICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#32CD32") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Household")


## office
estimate <- GTest[1,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

officeICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$office, function(x) arrange_ICCI(x$inferred.Tinfection)))
)


poffice <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = officeICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#808080") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Office work")


# Restaurant
estimate <- GTest[6,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

restaurantICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$restaurant, function(x) arrange_ICCI(x$inferred.Tinfection)))
)


pdine <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = restaurantICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#FF8C00") +
  xlab(" ") + ylab("Relative frequency") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Restaurant")


## Blue collar
estimate <- GTest[3,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

blueICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$manual_labour, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

pblue <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = blueICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#1E90FF") +
  xlab("") + ylab("Relative frequency") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Manual labour")


## GCLA
estimate <- GTest[2,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

retailICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$retail, function(x) arrange_ICCI(x$inferred.Tinfection)))
)
pGCLA <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = retailICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#008080") +
  xlab(" ") + ylab("") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Retail & leisure")
rm(GTest)

### Nosocomial
estimate <- GTest[7,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

nosoICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$noso, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

pnoso <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = nosoICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#800080") +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Nosocomial")
rm(GTest)

# indoor catering
estimate <- GTest[5,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

socialICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$social, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

pICFS <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = socialICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#FF1493") +
  xlab("Sampled infection intervals (days) from first to successive cases ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Close-social indoor")
rm(GTest)

### care homes
estimate <- GTest[8,]
shape.est = as.numeric(estimate[,2]^2/estimate[,5]^2)
rate.est = as.numeric(estimate[,2]/estimate[,5]^2)
w1 = as.numeric(estimate[,8])
w2 = as.numeric(estimate[,9])
w3 = as.numeric(estimate[,10])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[,2])^2/(2*estimate[,5]^2))
rate.sec = as.numeric(2*estimate[,2]/(2*estimate[,5]^2))
shape.tert = as.numeric((3*estimate[,2])^2/(3*estimate[,5]^2))
rate.tert = as.numeric(3*estimate[,2]/(3*estimate[,5]^2))

careICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$carehome, function(x) arrange_ICCI(x$inferred.Tinfection)))
)


pcare <- ggplot() +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y)) +
  geom_histogram(data = careICCI,
                 aes(x = ICCI, y = after_stat(density)), binwidth = 1, alpha = 0.5, color = "white", fill = "#8B4513") +
  xlab("Sampled infection intervals (days) from first to successive cases") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Care home")
rm(GTest)

library(ggpubr)

pdistcomb <- ggarrange(phousehold, poffice, pdine, pblue, pGCLA,
                       pnoso, pICFS, pcare, nrow = 4, ncol = 2, labels = letters[1:8]) 

ggsave("FigS5_20241212.pdf", pdistcomb, width = 10, height = 7.5, units = "in")

