## according to S. Pei et al., Dist of infection to test (D_itt)
empDist <- read.csv("SP_etal_infection_to_test_extract data.csv", header = F)
head(empDist)
empDist$Days <- round(empDist$V1)
empDist$PseudoFreq <- round(empDist$V2 * 1000)
empdays <- c()
length(empDist$Days)
for(i in 1:19){
  empdays <- c(empdays, rep(empDist$Days[i], empDist$PseudoFreq[i]))
}

library(fitdistrplus)

fit_Ditt <- fitdist(empdays, "gamma")
fit_Ditt$estimate

bootfit <- bootdist(fit_Ditt, niter = 100)
quantile(bootfit$estim[,1]/bootfit$estim[,2], c(0.05, 0.975)) # 7.81 8.20
quantile(sqrt(bootfit$estim[,1]/bootfit$estim[,2]^2), c(0.05, 0.975)) # 3.17 3.46



Ditt_shape <- fit_Ditt$estimate[1] # 5.83
Ditt_rate <- fit_Ditt$estimate[2] # 0.73


Ditt_shape/Ditt_rate # 7.98
sqrt(Ditt_shape/Ditt_rate^2) # 3.30

