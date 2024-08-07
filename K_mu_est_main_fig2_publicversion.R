library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

###
source("tupper.R")
load("size_by_setting_0410.rda")
clustall <- read.csv("allclustersize.csv")

outputres <- function(modelres){
  data.frame(
    dispersion_k = modelres$a_n_b.est[1],
    meanclust = 1 + modelres$a_n_b.est[1] * modelres$a_n_b.est[2],
    Rc = modelres$a_n_b.est[1] * modelres$a_n_b.est[2]
  )
}

## household size

famsize <- size_by_setting$household$clustsize

res.fam <- optim_fit_model(init = c(0.5, 0.5), 
                           inputdata = list(q1 = 0.9, # suppose q1 = 0.9
                                            q2 = 0.9, # suppose q2 = 0.9
                                            clusters = famsize))
res.fam$a_n_b.est[1] # 1.081795
1 + res.fam$a_n_b.est[1] * res.fam$a_n_b.est[2] # mean clust size = 2.12, Cz = 1.12

conf.famclust <- conf_ellipse(inputdata = list(a.est = res.fam$a_n_b.est[1],
                                               b.est = res.fam$a_n_b.est[2],
                                               hessmat = res.fam$hessian_mat),
                              conflevel = 0.95,
                              npoint = 100)
range(conf.famclust$x.ellipsepoints) #  2.057871 2.188999
range(conf.famclust$y.ellipsepoints) # 0.9146044 1.2489865


## care home size
caresize <- size_by_setting$carehome$clustsize

res.care <- optim_fit_model(init = c(0.5, 0.5), 
                            inputdata = list(q1 = 0.9, # suppose q1 = 0.9
                                             q2 = 0.6, # suppose q2 = 0.6
                                             clusters = caresize))
res.care$a_n_b.est[1] # 0.08385441
1 + res.care$a_n_b.est[1] * res.care$a_n_b.est[2] # mean clust size = 6.68

conf.careclust <- conf_ellipse(inputdata = list(a.est = res.care$a_n_b.est[1],
                                                b.est = res.care$a_n_b.est[2],
                                                hessmat = res.care$hessian_mat),
                               conflevel = 0.95,
                               npoint = 100)
range(conf.careclust$x.ellipsepoints) # -, 19.13
range(conf.careclust$y.ellipsepoints) # 0.02 0.14

## restaurant

dinesize <- size_by_setting$restaurant$clustsize

res.dine <- optim_fit_model(init = c(0.5, 0.5), 
                            inputdata = list(q1 = 0.8, # suppose q1 = 0.8
                                             q2 = 0.5, # suppose q2 = 0.5
                                             clusters = dinesize))
outputres(res.dine)
#   dispersion_k meanclust  C_Z
#    0.4235245  2.451511 1.451511
conf.dineclust <- conf_ellipse(inputdata = list(a.est = res.dine$a_n_b.est[1],
                                                b.est = res.dine$a_n_b.est[2],
                                                hessmat = res.dine$hessian_mat),
                               conflevel = 0.95,
                               npoint = 100)
range(conf.dineclust$x.ellipsepoints) #  2.120791 2.782230
range(conf.dineclust$y.ellipsepoints) # 0.2817965 0.5652524

## nosocomial

nososize <- size_by_setting$nosocomial$clustsize

res.noso <- optim_fit_model(init = c(0.5, 0.5), 
                            inputdata = list(q1 = 0.9, 
                                             q2 = 0.7, 
                                             clusters = nososize))

outputres(res.noso)
#  dispersion_k meanclust  C_Z
#    0.2404996  2.194468 1.194468
conf.nosoclust <- conf_ellipse(inputdata = list(a.est = res.noso$a_n_b.est[1],
                                                b.est = res.noso$a_n_b.est[2],
                                                hessmat = res.noso$hessian_mat),
                               conflevel = 0.95,
                               npoint = 100)
range(conf.nosoclust$x.ellipsepoints) # 1.489503 2.899434
range(conf.nosoclust$y.ellipsepoints) # 0.07609519 0.40490405

### blue collar (in manuscript we call it manual labour)
bluesize <- size_by_setting$blue$clustsize
res.blue <- optim_fit_model(init = c(0.5, 0.5), 
                            inputdata = list(q1 = 0.9, # suppose q1 = 0.9
                                             q2 = 0.7, # suppose q2 = 0.7
                                             clusters = bluesize))
outputres(res.blue)
# dispersion_k meanclust  C_Z
#    0.1270421  2.104296 1.104296
conf.blueclust <- conf_ellipse(inputdata = list(a.est = res.blue$a_n_b.est[1],
                                                b.est = res.blue$a_n_b.est[2],
                                                hessmat = res.blue$hessian_mat),
                               conflevel = 0.95,
                               npoint = 100)
range(conf.blueclust$x.ellipsepoints) # 1.648333 2.560259
range(conf.blueclust$y.ellipsepoints) # 0.07739427 0.17669003

### office work
## in code I denote as white-collar
whitesize <- size_by_setting$office$clustsize

res.white <- optim_fit_model(init = c(0.5, 0.5), 
                             inputdata = list(q1 = 0.9, # suppose q1 = 0.9
                                              q2 = 0.7, # suppose q2 = 0.7
                                              clusters = whitesize))
outputres(res.white)
#  dispersion_k meanclust CZ
#    0.2253797  1.324696 0.324696
conf.whiteclust <- conf_ellipse(inputdata = list(a.est = res.white$a_n_b.est[1],
                                                 b.est = res.white$a_n_b.est[2],
                                                 hessmat = res.white$hessian_mat),
                                conflevel = 0.95,
                                npoint = 100)
range(conf.whiteclust$x.ellipsepoints) # 1.211014 1.438378
range(conf.whiteclust$y.ellipsepoints) # 0.08677858 0.36398072

### GCLA general consumption & leisure activities
## finally we re-termed it as retail and leisure
GCLAsize <- size_by_setting$retail$clustsize
res.GCLA <- optim_fit_model(init = c(0.5, 0.5), 
                            inputdata = list(q1 = 0.8, # suppose q1 = 0.9
                                             q2 = 0.5, # suppose q2 = 0.5
                                             clusters = GCLAsize))
outputres(res.GCLA)
#  dispersion_k meanclust   CZ
#   0.05200852  1.631335 0.6313354
conf.GCLAclust <- conf_ellipse(inputdata = list(a.est = res.GCLA$a_n_b.est[1],
                                                b.est = res.GCLA$a_n_b.est[2],
                                                hessmat = res.GCLA$hessian_mat),
                               conflevel = 0.95,
                               npoint = 100)
range(conf.GCLAclust$x.ellipsepoints) # 1.071919 2.190751
range(conf.GCLAclust$y.ellipsepoints) # 0.01017502 0.09384202

### ICFS (indoor-catering-fitness-social)
## finally we re-termed as close-indoor and social
ICFSsize <- size_by_setting$indoor$clustsize
res.ICFS <- optim_fit_model(init = c(0.5, 0.5), 
                            inputdata = list(q1 = 0.8, # suppose q1 = 0.8
                                             q2 = 0.4, # suppose q2 = 0.4
                                             clusters = ICFSsize))
outputres(res.ICFS)
# dispersion_k meanclust  CZ
#  0.07076718  7.067132 6.067132
conf.ICFSclust <- conf_ellipse(inputdata = list(a.est = res.ICFS$a_n_b.est[1],
                                                b.est = res.ICFS$a_n_b.est[2],
                                                hessmat = res.ICFS$hessian_mat),
                               conflevel = 0.95,
                               npoint = 100)
range(conf.ICFSclust$x.ellipsepoints) #  1.969202 12.165062
range(conf.ICFSclust$y.ellipsepoints) # 0.03067296 0.11086141


### all clusters
clustsize.all <- clustall$clustsize
res.all <- optim_fit_model(init = c(0.5, 0.5),
                           inputdata = list(q1 = 0.8,
                                            q2 = 0.7,
                                            clusters = clustsize.all))
outputres(res.all)

# dispersion_k meanclust    CZ
#   0.4970277  2.478496 1.478496
conf.allclust <- conf_ellipse(inputdata = list(a.est = res.all$a_n_b.est[1],
                                               b.est = res.all$a_n_b.est[2],
                                               hessmat = res.all$hessian_mat),
                              conflevel = 0.95,
                              npoint = 100)
range(conf.allclust$x.ellipsepoints) #  2.386314 2.570678
range(conf.allclust$y.ellipsepoints) #  0.4492078 0.5448476

########## visualization of ellipse

conf.df <- data.frame(
  xconfpoints = c(conf.allclust$x.ellipsepoints,
                  conf.famclust$x.ellipsepoints,
                  conf.dineclust$x.ellipsepoints,
                  conf.ICFSclust$x.ellipsepoints,
                  conf.GCLAclust$x.ellipsepoints,
                  conf.blueclust$x.ellipsepoints,
                  conf.whiteclust$x.ellipsepoints,
                  conf.careclust$x.ellipsepoints,
                  conf.nosoclust$x.ellipsepoints
  ),
  yconfpoints = c(conf.allclust$y.ellipsepoints,
                  conf.famclust$y.ellipsepoints,
                  conf.dineclust$y.ellipsepoints,
                  conf.ICFSclust$y.ellipsepoints,
                  conf.GCLAclust$y.ellipsepoints,
                  conf.blueclust$y.ellipsepoints,
                  conf.whiteclust$y.ellipsepoints,
                  conf.careclust$y.ellipsepoints,
                  conf.nosoclust$y.ellipsepoints
  ),
  cluster = c(rep("All cluster", length(conf.allclust$x.ellipsepoints)),
              rep("Households", length(conf.famclust$x.ellipsepoints)),
              rep("Restaurants", length(conf.blueclust$x.ellipsepoints)),
              rep("Close-social indoor", length(conf.ICFSclust$x.ellipsepoints)),
              rep("Retail & leisure", length(conf.GCLAclust$x.ellipsepoints)),
              rep("Manual labour", length(conf.blueclust$x.ellipsepoints)),
              rep("Office work", length(conf.whiteclust$x.ellipsepoints)),
              rep("Care homes", length(conf.careclust$x.ellipsepoints)),
              rep("Nosocomial", length(conf.nosoclust$x.ellipsepoints)))
  
)

conf.df$cluster <- factor(conf.df$cluster, levels = c("All cluster", 
                                                      "Households",
                                                      "Office work",
                                                      "Restaurants",
                                                      "Manual labour",
                                                      "Retail & leisure",
                                                      "Nosocomial",
                                                      "Close-social indoor",
                                                      "Care homes"))



pointdf <- data.frame(
  x = c(outputres(res.all)[2]%>% as.numeric(), 
        outputres(res.fam)[2]%>% as.numeric(), 
        outputres(res.dine)[2]%>% as.numeric(),
        outputres(res.ICFS)[2]%>% as.numeric(),
        outputres(res.GCLA)[2]%>% as.numeric(),
        outputres(res.blue)[2]%>% as.numeric(),
        outputres(res.white)[2]%>% as.numeric(),
        outputres(res.care)[2]%>% as.numeric(), 
        outputres(res.noso)[2]%>% as.numeric()),
  y = c(outputres(res.all)[1]%>% as.numeric(), 
        outputres(res.fam)[1]%>% as.numeric(),
        outputres(res.dine)[1]%>% as.numeric(),
        outputres(res.ICFS)[1]%>% as.numeric(),
        outputres(res.GCLA)[1]%>% as.numeric(),
        outputres(res.blue)[1]%>% as.numeric(),
        outputres(res.white)[1]%>% as.numeric(),
        outputres(res.care)[1]%>% as.numeric(), 
        outputres(res.noso)[1]%>% as.numeric()),
  cluster = c("All cluster", "Households","Restaurants", "Close-social indoor", "Retail & leisure",
              "Manual labour","Office work","Care homes","Nosocomial")
)

pointdf$cluster <- factor(pointdf$cluster, levels = c("All cluster", 
                                                      "Households",
                                                      "Office work",
                                                      "Restaurants",
                                                      "Manual labour",
                                                      "Retail & leisure",
                                                      "Nosocomial",
                                                      "Close-social indoor",
                                                      "Care homes"))

ellipse.1 <- ggplot() + 
  geom_path(data = conf.df , 
            aes(x = xconfpoints, y = yconfpoints, group = cluster, color = cluster)) +
  geom_point(data = pointdf , 
             aes(x = x, y = y, group = cluster, color = cluster)) + 
  scale_color_manual(values=c("#000000", 
                              "#32CD32", 
                              "#808080",
                              "#FF8C00",
                              "#1E90FF",
                              "#008080",
                              "#800080",
                              "#FF1493",
                              "#8B4513"))+
  theme_bw() + theme(legend.position = "right") +
  scale_x_log10(limits = c(1, 20), breaks = c(1, 1.5, 2, 2.5, 3, 7)) +
  # scale_y_continuous(limits = c(0, 1.5)) +
  scale_y_log10(breaks = c(0.01, 0.05, 0.10, 0.25, 0.5, 1.0), limits = c(0.01, 1.5)) +
  ylab("Dispersion parameter, k") + 
  xlab(expression(paste("Mean cluster size (1 + ", C[Z], ")")))
ellipse.1


#### 20-80


propresponsible <- function(R0, k, prop = 0.8) {
  qm1 <- qnbinom(1 - prop, k + 1, mu = R0 * (k + 1) / k)
  remq <- 1 - prop - pnbinom(qm1 - 1, k + 1, mu = R0 * (k + 1) / k)
  remx <- remq / dnbinom(qm1, k + 1, mu = R0 * (k + 1) / k)
  q <- qm1 + 1
  1 - pnbinom(q - 1, k, mu = R0) - dnbinom(q, k, mu = R0) * remx
}



frac_est_CI.2080 <- function(confinput, estinput){
  frac_all = numeric(201)
  for(i in 1:201){
    frac_all[i] <- propresponsible(R0 = confinput$x.ellipsepoints[i] - 1,
                                   k = confinput$y.ellipsepoints[i], prop = 0.8)
  }
  frac_CI = range(frac_all)
  frac_est = propresponsible(R0 = as.numeric(estinput[3]),
                             k = as.numeric(estinput[1]),
                             prop = 0.8)
  return(list(
    est_fraction = frac_est,
    CI_fraction = frac_CI
  ))
}
# 
propCI_allclust <- frac_est_CI.2080(confinput = conf.allclust,
                                    estinput = outputres(res.all))
propCI_allclust

propCI_fam <- frac_est_CI.2080(confinput = conf.famclust,
                               estinput = outputres(res.fam))
propCI_fam

propCI_dine <- frac_est_CI.2080(confinput = conf.dineclust,
                                estinput = outputres(res.dine))
propCI_dine

propCI_GCLA <- frac_est_CI.2080(confinput = conf.GCLAclust,
                                estinput = outputres(res.GCLA))
propCI_GCLA

propCI_ICFS <- frac_est_CI.2080(confinput = conf.ICFSclust,
                                estinput = outputres(res.ICFS))
propCI_ICFS

propCI_blue <- frac_est_CI.2080(confinput = conf.blueclust,
                                estinput = outputres(res.blue))
propCI_blue

propCI_white <- frac_est_CI.2080(confinput = conf.whiteclust,
                                 estinput = outputres(res.white))
propCI_white

propCI_noso <- frac_est_CI.2080(confinput = conf.nosoclust,
                                estinput = outputres(res.noso))
propCI_noso

### for carehome clusters, remove the points where CZ is smaller than 0
ind.exc <- which(conf.careclust$x.ellipsepoints <= 0)
### note the function is generalized to CZ, just in previous version we termed as R0
### but CZ is equal or larger than R0

propresponsible(R0 = as.numeric(outputres(res.care)[3]),
                k = as.numeric(outputres(res.care)[1]),
                prop = 0.8) # 0.08232631

201 - length(ind.exc) # 138
frac_careclust <- numeric(138)
for(i in 1:138){
  xpoint <- conf.careclust$x.ellipsepoints[i]
  ypoint <- conf.careclust$y.ellipsepoints[i]
  if(xpoint < 1){xpoint = 1.001}
  frac_careclust[i] <- propresponsible(R0 = xpoint - 1,
                                       k = ypoint, prop = 0.8)
}
range(frac_careclust)


prop_20_80_list <- list(
  allclust = propCI_allclust,
  famclust = propCI_fam,
  careclust = range(frac_careclust),
  dineclust = propCI_dine,
  nosoclust = propCI_noso,
  blueclust = propCI_blue,
  whiteclust = propCI_white,
  indoorclust = propCI_ICFS,
  GCLAclust = propCI_GCLA
)

# save(prop_20_80_list, file = "prop_20_80_list.rda")


Rvect <- seq(0.01, 10, 0.01)
kvect <- seq(0.01, 1.5, 0.001)

p80 <- matrix(
  nrow = length(kvect),
  ncol = length(Rvect),
  byrow = T
)

# load("prop_20_80_list.rda")

propresponsible <- function(R0, k, prop = 0.8) {
  qm1 <- qnbinom(1 - prop, k + 1, mu = R0 * (k + 1) / k)
  remq <- 1 - prop - pnbinom(qm1 - 1, k + 1, mu = R0 * (k + 1) / k)
  remx <- remq / dnbinom(qm1, k + 1, mu = R0 * (k + 1) / k)
  q <- qm1 + 1
  1 - pnbinom(q - 1, k, mu = R0) - dnbinom(q, k, mu = R0) * remx
}


for (jj in 1:length(Rvect)) {
  R <- Rvect[jj]
  for (kk in 1:length(kvect)) {
    k <- kvect[kk]
    tmp <- propresponsible(R, k, prop = 0.8)
    p80[kk, jj] <- tmp[[1]]
  }
}


mycontour <- data.frame(
  x = rep(Rvect, each = length(kvect)) + 1,
  y = rep(kvect, length(Rvect)),
  p80_out = c(p80)
)

pcontour <- ggplot(mapping = aes(x, y)) + 
  geom_contour_filled(data = mycontour, aes(z = p80_out, fill = after_stat(level)), 
                      breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.5),
                      color = "black", linewidth = 0.25) + 
  scale_fill_brewer(palette = "OrRd", direction = -1) +
  scale_y_log10(breaks = c(0.01, 0.05, 0.10, 0.25, 0.5, 1.0), limits = c(0.01, 1.5)) +
  #  scale_x_continuous(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 6.5, 7.0, 7.5), limits = c(1, 8)) + 
  scale_x_log10(limits = c(1, 20), breaks = c(1, 1.5, 2, 2.5, 3, 7)) +
  geom_point(data = data.frame(x = pointdf$x[1], y = pointdf$y[1]), size = 1.5) +
  annotate("text", x = pointdf$x[1] + 0.65, y = pointdf$y[1], size = 3.5, label = "All clusters") +
  
  geom_point(data = data.frame(x = pointdf$x[2], y = pointdf$y[2]), size = 1.5) +
  annotate("text", x = pointdf$x[2] + 0.7, y = pointdf$y[2], size = 3.5, label = "Households") +
  
  geom_point(data = data.frame(x = pointdf$x[3], y = pointdf$y[3]), size = 1.5) +
  annotate("text", x = pointdf$x[3] + 0.65, y = pointdf$y[3] - 0.01, size = 3.5, label = "Restaurants") +
  
  geom_point(data = data.frame(x = pointdf$x[4], y = pointdf$y[4]), size = 1.5) +
  annotate("text", x = pointdf$x[4] , y = pointdf$y[4] - 0.01, size = 3.5, label = "Close-social indoor") +
  
  geom_point(data = data.frame(x = pointdf$x[5], y = pointdf$y[5]), size = 1.5) +
  annotate("text", x = pointdf$x[5] + 0.35, y = pointdf$y[5] - 0.01, size = 3.5, label = "Retail & leisure") +
  
  geom_point(data = data.frame(x = pointdf$x[6], y = pointdf$y[6]), size = 1.5) +
  annotate("text", x = pointdf$x[6] + 0.75, y = pointdf$y[6], size = 3.5, label = "Manual labour") +
  
  geom_point(data = data.frame(x = pointdf$x[7], y = pointdf$y[7]), size = 1.5) +
  annotate("text", x = pointdf$x[7] + 0.25, y = pointdf$y[7] - 0.025, size = 3.5, label = "Office work") +
  
  geom_point(data = data.frame(x = pointdf$x[8], y = pointdf$y[8]), size = 1.5) +
  annotate("text", x = pointdf$x[8] + 1.75, y = pointdf$y[8], size = 3.5, label = "Care homes") +
  
  geom_point(data = data.frame(x = pointdf$x[9], y = pointdf$y[9]), size = 1.5) +
  annotate("text", x = pointdf$x[9] + 0.6, y = pointdf$y[9], size = 3.5, label = "Nosocomial") +
  
  labs(x = expression(paste("Mean cluster size (1 + ", C[Z], ")")), 
       y = "Dispersion parameter, k", 
       fill = "Prop. 80%") +
  
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        legend.position = c(0.9, 0.5),
        panel.grid = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank()) 


pfig2_main <- ggarrange(ellipse.1, pcontour, nrow = 1, ncol = 2, align = "h", labels = c("a", "b"))
pfig2_main


# ggsave("Fig2_main.pdf", pfig2_main, width = 13, height = 6, units = "in", dpi = 900)


























