library(readxl)
library(scales)
library(ggpubr)

### from data/processed
estGT <- read_excel("est_GI_combine.xlsx")

estGT$cluster <- c("Households", "Care homes", "Restaurants", "Nosocomial",
                   "Manual labour", "Office work", "Close-social indoor", "Retail & leisure")
estGT$cluster <- factor(estGT$cluster, levels = c("Households", "Office work", "Restaurants", "Manual labour",
                                                  "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes"))

estGT$X <- c(1, 8, 3, 6, 4, 2, 7, 5)/2


plotmu <- ggplot() + geom_point(data = estGT, aes(x = X, y = mu, color = cluster), size = 2) +
  geom_errorbar(data = estGT, aes(x = X, ymin = muLB, ymax = muUB, color = cluster), size = 1, width = 0) +
  scale_color_manual(values = c("#32CD32", 
                                "#808080", "#FF8C00",
                                "#1E90FF", "#008080",
                                "#800080", "#FF1493", 
                                "#8B4513")) +
  scale_x_continuous(breaks = seq(0.5, 4, 0.5),
                     labels = c("Households", "Office work", "Restaurants", "Manual labour",
                                "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes")) +
  scale_y_continuous(breaks = seq(2.5, 10, 2.5),
                     labels = label_number(accuracy = 0.1)) +
  geom_hline(yintercept = 5.6, linetype = "dashed") +
  xlab(" ") + ylab("Days") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "none",
        axis.ticks.x.bottom = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ggtitle("Mean generation interval") 

estGT$cluster

plotsd <- ggplot() + geom_point(data = estGT, aes(x = X, y = sd, color = cluster), size = 2) +
  geom_errorbar(data = estGT, aes(x = X, ymin = sdLB, ymax = sdUB, color = cluster), size = 1, width = 0) +
  scale_color_manual(values = c("#32CD32", 
                                "#808080", "#FF8C00",
                                "#1E90FF", "#008080",
                                "#800080", "#FF1493", 
                                "#8B4513")) +
  scale_x_continuous(breaks = seq(0.5, 4, 0.5),
                     labels = c("Households", "Office work", "Restaurants", "Manual labour",
                                "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes")) +
  scale_y_continuous(breaks = seq(0.5, 4.0, 1.5)) +
  xlab(" ") + ylab("Days") + 
  geom_hline(yintercept = 1.9, linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ggtitle("Standard deviation of generation interval")



estGT$shape <- estGT$mu^2/estGT$sd^2
estGT$rate <- estGT$mu/estGT$sd^2


plotpdf <- ggplot() + 
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[1],
                                         rate = estGT$rate[1])},
                aes(color = paste(2, "Households")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[2],
                                         rate = estGT$rate[2])},
                aes(color = paste(9, "Care homes")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[3],
                                         rate = estGT$rate[3])},
                aes(color = paste(4, "Restaurants")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[4],
                                         rate = estGT$rate[4])},
                aes(color = paste(7, "Nosocomial")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[5],
                                         rate = estGT$rate[5])},
                aes(color = paste(5, "Manual labour")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[6],
                                         rate = estGT$rate[6])},
                aes(color = paste(3, "Office work")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[7],
                                         rate = estGT$rate[7])},
                aes(color = paste(8, "Close-social indoor")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = estGT$shape[8],
                                         rate = estGT$rate[8])},
                aes(color = paste(6, "Retail & leisure")), size = 1) +
  geom_function(fun = function(x){dgamma(x, shape = 5.6^2/1.9^2,
                                         rate = 5.6/1.9^2)},
                aes(color = paste(1, "All pooled")), size = 1, linetype = "dashed")+
  scale_color_manual(
    name = c("Cluster setting"),
    values =  c("black", "#32CD32", 
                "#808080", "#FF8C00",
                "#1E90FF", "#008080",
                "#800080", "#FF1493", 
                "#8B4513"),
    labels = c("All pooled","Households", "Office work", "Restaurants", "Manual labour",
               "Retail & leisure", "Nosocomial", "Close-social indoor", "Care homes")) +
  scale_x_continuous(limits = c(0, 20)) +
  xlab("Days") + ylab("Probability") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        legend.position = c(0.825, 0.625),
        legend.title = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  ggtitle("Probability density of inferred generation interval distributions")

plotpdf

fig3 <- ggarrange(plotmu, plotsd, plotpdf, nrow = 3, ncol = 1, labels = letters[1:3], align = "v")


# ggsave("Fig3_main.pdf", fig3, width = 12, height = 8.5, units = "in", dpi = 900)
