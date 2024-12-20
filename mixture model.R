#adapted from WEB APPENDIX 2 from Vink, Bootsma & Wallinga 2014
library(fdrtool)


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

# # 
# pfgd <- function(x,a,b){
#   pt <- function(xi,a,b){
#     int <- integrate(f=dfgd, lower=0, upper=xi, a=a, b=b)
#     val <- int$value
#     return(val)
#   }
#   Ft <- sapply(1:length(x), function(i) pt(x[i], a, b))
#   Ft[which(Ft>1)] <- 1
#   return(Ft)
# }



# E-step
GI_mix_est <- function(data, N, startmu = mean(data), startsig = sd(data)) {
  
  # transfer to shape and rate
  a.pri = (startmu^2)/(startsig^2) 
  b.pri = startmu/(startsig^2)
  
  
  #function for weighted variance, to be used in E-step
  #weighted variance
  weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
      w <- w[i <- !is.na(x)]
      x <- x[i]
    }
    sum.w <- sum(w)
    (sum.w*(sum(w*(x-weighted.mean(x,w))^2))) / (sum.w^2 - sum(w^2))
  }
  
  #mixture with 4 components
  #we split the gamma (and folded gamma difference) distribution for the PS, PT and PQ route into two parts
  #component 1: CP route
  #component 2: PS route
  #component 3: PT route
  #component 4: PQ route
  
  # convolution of the triangular distribution with the mixture component density
  # contiunous case. This is bc we have only day of infection, not exact time. 
  
  # coprimary transmission: both infected by the same primary case. |x2-x1| -> folded gamma difference dist
  
  f10 <- function(x, a.pri, b.pri)    (2-2*x)*dfgd(x, a = a.pri, b = b.pri)
  f1lower <- function(x,r,a.pri, b.pri)    	(x-r+1)*dfgd(x, a = a.pri, b = b.pri)
  f1upper <- function(x,r,a.pri, b.pri)    	(r+1-x)*dfgd(x, a = a.pri, b = b.pri)
  
  
  # primary - secondary. here it follows a gamma distribution 
  f20<-function(x,a.pri, b.pri)      	(2-2*x)*dgamma(x, shape = a.pri, rate = b.pri)
  f2lower<-function(x,r,a.pri, b.pri) 	(x-r+1)*dgamma(x, shape = a.pri, rate = b.pri)
  f2upper<-function(x,r,a.pri, b.pri) 	(r+1-x)*dgamma(x, shape = a.pri, rate = b.pri)
  
  # Primary-tertiary transmission, gamma dist with mean = 2*mu; variance = 2*sigma^2
  # transfer to shape and rate
  a.sec = ((2*startmu)^2)/(2*startsig^2) 
  b.sec = 2*startmu/(2*startsig^2)
  
  
  f30<-function(x, a.sec, b.sec)       	(2-2*x)*dgamma(x, shape = a.sec, rate = b.sec)
  f3lower<-function(x,r,a.sec, b.sec) 	(x-r+1)*dgamma(x, shape = a.sec, rate = b.sec)
  f3upper<-function(x,r,a.sec, b.sec) 	(r+1-x)*dgamma(x, shape = a.sec, rate = b.sec)
  
  # Primary-quadternary transmission, gamma dist with mean = 3*mu; variance = 3*sigma^2
  # transfer to shape and rate
  a.tert = ((3*startmu)^2)/(3*startsig^2) 
  b.tert = 3*startmu/(3*startsig^2)
  
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

    
  j <- length(data) 
  
  #EM algorithm: 
  tau1 <- numeric(j)    #mixture weights for (coprimary, coprimary) pairs
  tau2 <- numeric(j)    #mixture weights for (primary, secondary) pairs
  tau3 <- numeric(j)    #mixture weights for (primary, tertiary) pairs 
  tau4 <- numeric(j)    #mixture weights for (primary, quaternery) pairs
  
  #number of iterations
  
  #calculate the absolute probability of interval belonging to a component
  for(k in 1:N){
    for(l in 1:j){
      if(data[l]==0){
        d1 <- p10(data[l], a.pri, b.pri)[[1]]
        d2 <- p20(data[l], a.pri, b.pri)[[1]]
        d3 <- p30(data[l], a.sec, b.sec)[[1]]
        d4 <- p40(data[l], a.tert, b.tert)[[1]]
      }
      else{
        d1<-p1lower(data[l], a.pri, b.pri)[[1]] + p1upper(data[l], a.pri, b.pri)[[1]]
        d2<-p2lower(data[l], a.pri, b.pri)[[1]] + p2upper(data[l], a.pri, b.pri)[[1]]
        d3<-p3lower(data[l], a.sec, b.sec)[[1]] + p3upper(data[l], a.sec, b.sec)[[1]]
        d4<-p4lower(data[l], a.tert, b.tert)[[1]] + p4upper(data[l], a.tert, b.tert)[[1]]
      }
      
      #then calculate the relative probability of a data point to belong to one of the components
      dummy<-d1+d2+d3+d4
      tau1[l]<-d1/dummy
      tau2[l]<-d2/dummy
      tau3[l]<-d3/dummy
      tau4[l]<-d4/dummy
    }
    
    #now calculate the weights for each of the components
    w1 <- sum(tau1)/j
    w2 <- sum(tau2)/j
    w3 <- sum(tau3)/j
    w4 <- sum(tau4)/j
    
    # M-step
    # estimates for the mean and standard deviation of the primary-secondary infection component can #be calculated directly
    mu <- weighted.mean(data, tau2)
    sigma <- sqrt(weighted.var(data, tau2))
    print(c(mu, sigma, w1, w2, w3))
    # update parameters in the loop
    a.pri = (mu^2)/(sigma^2) 
    b.pri = mu/(sigma^2)
    a.sec = ((2*mu)^2)/(2*sigma^2) 
    b.sec = 2*mu/(2*sigma^2)
    a.tert = ((3*mu)^2)/(3*sigma^2) 
    b.tert = 3*mu/(3*sigma^2)
  }
  return(c(mu, sigma, w1, w2, w3))
}
