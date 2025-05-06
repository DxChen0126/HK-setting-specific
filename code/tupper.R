# % Y is number of new infections is Poisson(gamma) where gamma is Gamma(a=k,b=theta)
# % X is binom(Y+1,q1)
# % then it is censored with probabilty 1-(1-q2)^X
# % output is vector where kth entry is Prob(X=k), k >= 1 
# % a is shape for gamma dist, b is scale, q is ascertainment probability
# function prob=generate_prob_vector(a,b,q1,q2,maxk)

generate_prob_vector <- function(a, b, q1, q2, maxk){
# % how big we have to go for Y depends on how big we go for X
# if q1<0
# error('q1 cannot be less than zero')
# end
# if q2<0
# error('q2 cannot be zero')
# end
# 
 maxi = round(maxk/q1 * 2)
# % Y is negative binomial with params r,p
 r = a
 p = b/(1 + b)
# % convert to q for matlab (actually I don't get this line, but R has the same order of para pdf as in matlab)
 q = 1 - p
# Yprob=nbinpdf((0:maxi),r,q);
 Yprob = dnbinom(0:maxi, r, q)
# % Yprob(i)=Prob(Y+1 = i)
# prob=zeros(maxk,1);
 prob = numeric(maxk) 
# % k is the value of X
# for k=0:maxk
 for(k in 0:maxk){
# % range of possible Y that could give X=k
   ivect = max(k-1, 0):maxi
# ivect=(max([k-1 0]):maxi);
# % prob X=k given Y=i
   probXgivenY = dbinom(k, ivect+1, q1)
# probXgivenY=binopdf(k,ivect+1,q1);
# % prob that X is k
   prob[k+1] = sum(probXgivenY * Yprob[ivect+1])
# prob(k+1)=sum( probXgivenY.*Yprob(ivect+1) );
# end
}
# % censor with probability 1-(1-q2)^X
# prob=prob.*( 1 - (1-q2).^(0:maxk)' );
 prob = prob * (1 - (1-q2) ^ t(0:maxk))
# % omit the first term (meaning zero observed infections) and normalize to get probability
# prob=prob(2:end)/sum(prob(2:end));
# I think in our situation is that we did not observe cluster size of 1
# for specific events too?
 prob = prob[-1]/sum(prob[-1]) 
# 
# end
}


# negative log likelihood 
negLogLik <- function(pars, data){
  # pars means the parameter a and b to be estimated
  # data is a list that contains choices of q1, q2, and cluster vector
  # cluster vector means a vactor of cluster distributions
  a = pars[1]
  b = pars[2]
  q1 = data[[1]]
  q2 = data[[2]]
  clusters = data[[3]]
  max_i = max(clusters)
  prob = generate_prob_vector(a, b, q1, q2, max_i)
  return(-sum(log(prob[clusters]))) # [clusters+1] tp avoid prob 0
}



optim_fit_model <- function(init, inputdata){
  # initial values for parameters a, b
  # a vector for example c(0.5, 0.5)
  res = optim(par = init, 
               fn = negLogLik,
               method = "Nelder-Mead",
               hessian = T,
               data = inputdata)
  pars.est = res$par
  negsumllk = res$value
  hessmat = res$hessian
  return(list(
    a_n_b.est = pars.est,
    negsumllk = negsumllk,
    hessian_mat = res$hessian
  ))
}


## suppose the output model fit is modeloutput
## modeloutput <- optim_fit_model(iniit = c(0.5, 0.5), inputdata = ...)
# a.est <- modeloutput$a_n_b.est[1]
# b.est <- modeloutput$a_n_b.est[2]
# meanClust <- 1 + a.est * b.est
# kDispersion <- a.est

# testca <- read.csv("testdata.csv")
# testres <- optim_fit_model(init = c(0.5, 0.5), 
#                            inputdata = list(q1 = 0.75,
#                                             q2 = 1,
#                                             clusters = testca$X1))


# it seems to work, at least can replicate Tupper's results.
# the problem is in our data we lack cluster size = 1, thus making the estimate of R uncertain
# maybe look into those single case profiles, could assign some single cases to
# other events according to the case work and building profile (if available)

### confidence ellipse

conf_ellipse <- function(inputdata, conflevel, npoint){
  #inputdata is a list that correspond to testres results
  a.est = inputdata[[1]]
  b.est = inputdata[[2]]
  hessmat = inputdata[[3]] # hessian matrix
  # hessmat = testres$hessian_mat
  # a.est = a.test
  # b.est = b.test
  Jac = matrix(c(b.est, a.est, 1, 0), nrow = 2, ncol = 2, byrow = T)
  # covariance matrix
  Cmat = Jac %*% solve(hessmat) %*% t(Jac)
  # to make it consistent with matlab
  eigenval = eigen(Cmat)$values
  eigenvalmat = matrix(c(eigenval[2], 0, 0, eigenval[1]), nrow = 2, ncol = 2, byrow = T)
  eigenvec = eigen(Cmat)$vectors
  eigenvecmat = matrix(c(eigenvec[1,2], eigenvec[1,1], eigenvec[2,2], eigenvec[2,1]), 
                       nrow = 2, ncol = 2, byrow = T)
  ## generate x and y points that define an ellipse, given covariance matrix
  # angles around a circle
  # default npoint is 100
  p = seq(from = 0, to = 2*pi, by = pi/npoint)
  pimat = matrix(c(cos(p), sin(p)), ncol = 2, byrow = F)
  xy = pimat %*% sqrt(eigenvalmat) %*% t(eigenvecmat)
  xpoints = xy[,1]
  ypoints = xy[,2]
  
  ## confidence interval
  # Compute quantile for the desired percentile
  # conflevel = 0.95
  k = sqrt(qchisq(conflevel, df = 2))
  xpoints.conf = 1 + a.est * b.est + xpoints * k
  ypoints.conf = a.est + ypoints * k
  return(list(
    x.ellipsepoints = xpoints.conf,
    y.ellipsepoints = ypoints.conf,
    x.originpoints = xpoints,
    y.originpoints = ypoints
  ))
}

# bravo! get the confidence ellipse










