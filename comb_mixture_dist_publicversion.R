load("./sample infection record/present_GI_100times_sample_20231012.rda")

# titlenames <- c("Households", "Care homes", "Restaurants", "Nosocomial",
#                 "Blue-collar work", "Office work", "Indoor catering, fitness & social",
#                 "Retail & leisure activities")
# 
# pooltb <- sampleGT_record[[1]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
# pooltb
# colnames(pooltb) <- c("ICC interval" , "Frequency")



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

GIest <- read.csv("pooled_GI_est_main.csv")
colnames(GIest)

estimate <- GIest[1, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_household <- sampleGT_record[[1]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_household) <- c("ICC", "Frequency")
ICC_household$relative_fre <- ICC_household$Frequency/sum(ICC_household$Frequency)
ICC_household$ICC <- seq(0, 30, 1)

phousehold <- ggplot() +
  geom_rect(data = ICC_household,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5, color = "black", fill = "#32CD32", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  
  xlab(" ") + ylab("Relative frequency") +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Households")

phousehold

## office
estimate <- GIest[6, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_office <- sampleGT_record[[6]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_office) <- c("ICC", "Frequency")
ICC_office$relative_fre <- ICC_office$Frequency/sum(ICC_office$Frequency)
ICC_office$ICC <- c(seq(0, 26, 1), 28)

poffice <- ggplot() +
  geom_rect(data = ICC_office,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5,  color = "black", fill = "#808080", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  xlab(" ") + ylab("Relative frequency") +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Office work")
poffice

# dine
estimate <- GIest[3, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_dine <- sampleGT_record[[3]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_dine) <- c("ICC", "Frequency")
ICC_dine$relative_fre <- ICC_dine$Frequency/sum(ICC_dine$Frequency)
ICC_dine$ICC <- seq(0, 30, 1)

pdine <- ggplot() +
  geom_rect(data = ICC_dine,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5,  color = "black", fill = "#FF8C00", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  xlab(" ") + ylab("Relative frequency") +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Restaurants")
pdine

## blue collar
estimate <- GIest[5, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_blue <- sampleGT_record[[5]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_blue) <- c("ICC", "Frequency")
ICC_blue$relative_fre <- ICC_blue$Frequency/sum(ICC_blue$Frequency)
ICC_blue$ICC <- seq(0, 30, 1)

pblue <- ggplot() +
  geom_rect(data = ICC_blue,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5, color = "black", fill = "#1E90FF", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  xlab(" ") + ylab("Relative frequency") +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Blue-collar")
pblue

### GCLA
estimate <- GIest[8, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_GCLA <- sampleGT_record[[8]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_GCLA) <- c("ICC", "Frequency")
ICC_GCLA$relative_fre <- ICC_GCLA$Frequency/sum(ICC_GCLA$Frequency)
ICC_GCLA$ICC <- seq(0, 30, 1)

pGCLA <- ggplot() +
  geom_rect(data = ICC_GCLA,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5, color = "black", fill = "#008080", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  xlab(" ") + ylab("Relative frequency") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Retail & leisure")
pGCLA

## nosocomial
estimate <- GIest[4, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_noso <- sampleGT_record[[4]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_noso) <- c("ICC", "Frequency")
ICC_noso$relative_fre <- ICC_noso$Frequency/sum(ICC_noso$Frequency)
ICC_noso$ICC <- seq(0, 30, 1)

pnoso <- ggplot() +
  geom_rect(data = ICC_noso,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5,  color = "black", fill = "#800080", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  xlab(" ") + ylab("Relative frequency") +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) + 
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Nosocomial")
pnoso

### indoor
estimate <- GIest[7, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_ICFS <- sampleGT_record[[7]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_ICFS) <- c("ICC", "Frequency")
ICC_ICFS$relative_fre <- ICC_ICFS$Frequency/sum(ICC_ICFS$Frequency)
ICC_ICFS$ICC <- seq(0, 30, 1)

pICFS <- ggplot() +
  geom_rect(data = ICC_ICFS,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5, color = "black", fill = "#FF1493", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  xlab("Sampled infection intervals (days) from first to successive cases") + ylab("Relative frequency") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Close-social indoor")
pICFS

### care home
estimate <- GIest[2, -1]
shape.est = as.numeric(estimate[1]^2/estimate[4]^2) 
rate.est = as.numeric(estimate[1]/estimate[4]^2)
w1 = as.numeric(estimate[7])
w2 = as.numeric(estimate[10])
w3 = as.numeric(estimate[13])
w4 = 1 - w1 - w2 - w3
shape.sec = as.numeric((2*estimate[1])^2/(2*estimate[4]^2))
rate.sec = as.numeric(2*estimate[1]/(2*estimate[4]^2))
shape.tert = as.numeric((3*estimate[1])^2/(3*estimate[4]^2))
rate.tert = as.numeric(3*estimate[1]/(3*estimate[4]^2))

ICC_care <- sampleGT_record[[2]][seq(1, 200, 2)] %>% unlist() %>% table() %>% as.data.frame()
colnames(ICC_care) <- c("ICC", "Frequency")
ICC_care$relative_fre <- ICC_care$Frequency/sum(ICC_care$Frequency)
ICC_care$ICC <- seq(0, 30, 1)

pcare <- ggplot() +
  geom_rect(data = ICC_care,
            aes(xmin = ICC - 0.5, xmax = ICC + 0.5, 
                ymin = 0, ymax = relative_fre), alpha = 0.5, color = "black", fill = "#8B4513", lwd = 0.5) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dfgd(seq(0, 30, 0.01), a = shape.est, b = rate.est)*w1),
            aes(x = x, y = y), lwd = 0.75) +
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.est, rate = rate.est)*w2),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.sec, rate = rate.sec)*w3),
            aes(x = x, y = y), lwd = 0.75) + 
  geom_line(data = data.frame(x = seq(0, 30, 0.01), 
                              y = dgamma(seq(0, 30, 0.01), shape = shape.tert, rate = rate.tert)*w4),
            aes(x = x, y = y), lwd = 0.75) +
  scale_y_continuous(limits = c(0, 0.107), breaks = seq(0, 0.1, 0.025)) +
  xlab("Sampled infection intervals (days) from first to successive cases") + ylab("Relative frequency") +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y.left = element_text(size = 8)) + ggtitle("Care homes")
pcare

pdistcomb <- ggarrange(phousehold, poffice, pdine, pblue, pGCLA,
                       pnoso, pICFS, pcare,
                       nrow = 4, ncol = 2, labels = letters[1:8]) 

# ggsave("./event cluster analysis/manuscript/FigS5_20240116.pdf", pdistcomb, width = 12, height = 9, units = "in")


