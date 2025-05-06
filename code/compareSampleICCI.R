load("sample_times_list.rda")
load("sample_times_list_alt.rda")


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

householdICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$household, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

householdICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$household, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

phousehold <- ggplot() +
  geom_histogram(data = householdICCI,
               aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = householdICCI_alt,
               aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Households")


officeICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$office, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

officeICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$office, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

poffice <- ggplot() +
  geom_histogram(data = officeICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = officeICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Office work")


restaurantICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$restaurant, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

restaurantICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$restaurant, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

prestaurant <- ggplot() +
  geom_histogram(data = restaurantICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = restaurantICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab(" ") + ylab("Relative frequency") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Restaurants")

blueICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$manual_labour, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

blueICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$manual_labour, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

pblue <- ggplot() +
  geom_histogram(data = blueICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = blueICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab(" ") + ylab("Relative frequency") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Manual labour")


retailICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$retail, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

retailICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$retail, function(x) arrange_ICCI(x$inferred.Tinfection)))
)


pretail <- ggplot() +
  geom_histogram(data = retailICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = retailICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Retail & leisure")

nosoICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$noso, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

nosoICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list$noso, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

pnoso <- ggplot() +
  geom_histogram(data = nosoICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = nosoICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab(" ") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Nosocomial")

socialICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$social, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

socialICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$social, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

psocial <- ggplot() +
  geom_histogram(data = socialICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = socialICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab("Sampled infection intervals (days) from first to successive cases") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Close-social indoor")


careICCI <- data.frame(
  ICCI = unlist(lapply(sample_times_list$carehome, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

careICCI_alt <- data.frame(
  ICCI = unlist(lapply(sample_times_list_alt$carehome, function(x) arrange_ICCI(x$inferred.Tinfection)))
)

pcare <- ggplot() +
  geom_histogram(data = careICCI,
                 aes(x = ICCI, y = after_stat(density)), fill = "blue", alpha = 0.2, binwidth = 1) +
  geom_histogram(data = careICCI_alt,
                 aes(x = ICCI, y = after_stat(density)), fill = "red", alpha = 0.2, binwidth = 1) +
  xlab("Sampled infection intervals (days) from first to successive cases") + ylab(" ") +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_bw() + theme(legend.title = element_blank(),
                     panel.grid.minor = element_blank()) + ggtitle("Care homes")

pICCcomb <- ggarrange(phousehold, poffice, prestaurant, pblue, pretail,
                       pnoso, psocial, pcare, nrow = 4, ncol = 2, labels = letters[1:8]) 

ggsave("FigS3_20241212.pdf", pICCcomb, width = 12, height = 7.5, units = "in")





