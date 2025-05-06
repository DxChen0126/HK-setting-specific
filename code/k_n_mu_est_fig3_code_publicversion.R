library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

source("tupper.R")

load("all_cluster_comb_df.rda")
load("clustersize_list.rda")

outputres <- function(modelres){
  data.frame(
    dispersion_k = modelres$a_n_b.est[1],
    meanclust = 1 + modelres$a_n_b.est[1] * modelres$a_n_b.est[2],
    Rc = modelres$a_n_b.est[1] * modelres$a_n_b.est[2]
  )
}

### note in code notation we did not distinguish Rc (setting-specific R) and 
### Cz (mean new infection per cluster)

### all cluster
allclustsize <- allclustcomb.df %>% filter(category == "All clusters")
allclustsize <- allclustsize$size

allclust.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.8, 
                                                   q2 = 1, 
                                                   clusters = allclustsize))

outputres(allclust.mod1)
# dispersion_k meanclust       Rc
# 0.6339083  2.765101 1.765101
conf.allclust.mod1 <- conf_ellipse(inputdata = list(a.est = allclust.mod1$a_n_b.est[1],
                                                          b.est = allclust.mod1$a_n_b.est[2],
                                                          hessmat = allclust.mod1$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.allclust.mod1$x.ellipsepoints) # 2.660580 2.869623
range(conf.allclust.mod1$y.ellipsepoints) # 0.5735350 0.6942817

# cluster ascertainment, q1 = 1, give q2 a value
allclust.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 1, 
                                                   q2 = 0.8, 
                                                   clusters = allclustsize))

outputres(allclust.mod2)
# dispersion_k meanclust        Rc
# 0.5603333  2.294081 1.294081
conf.allclust.mod2 <- conf_ellipse(inputdata = list(a.est = allclust.mod2$a_n_b.est[1],
                                                          b.est = allclust.mod2$a_n_b.est[2],
                                                          hessmat = allclust.mod2$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.allclust.mod2$x.ellipsepoints) #  2.214803 2.373359
range(conf.allclust.mod2$y.ellipsepoints) # 0.5092295 0.6114371

# double ascertainment
allclust.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.8, 
                                                   q2 = 0.8, 
                                                   clusters = allclustsize))

outputres(allclust.mod3)
# dispersion_k meanclust     Rc
# 1     0.5412194  2.567919 1.567919
conf.allclust.mod3 <- conf_ellipse(inputdata = list(a.est = allclust.mod3$a_n_b.est[1],
                                                          b.est = allclust.mod3$a_n_b.est[2],
                                                          hessmat = allclust.mod3$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.allclust.mod3$x.ellipsepoints) #   2.471651 2.664187
range(conf.allclust.mod3$y.ellipsepoints) # 0.4894451 0.5929936

##################################
######### household setting
head(allclustcomb.df)
unique(allclustcomb.df$category)
householdsize <- allclustcomb.df %>% filter(category == "Households") 
householdsize <- householdsize$size
# individual ascertainment, q2 = 1, give q1 a value
household.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 0.9, 
                                                q2 = 1, 
                                                clusters = householdsize))

outputres(household.mod1)
# dispersion_k meanclust       Rc
# 1.257068  2.169965 1.169965
conf.householdclust.mod1 <- conf_ellipse(inputdata = list(a.est = household.mod1$a_n_b.est[1],
                                                    b.est = household.mod1$a_n_b.est[2],
                                                    hessmat = household.mod1$hessian_mat),
                                   conflevel = 0.95,
                                   npoint = 100)

range(conf.householdclust.mod1$x.ellipsepoints) #  2.103883 2.236047
range(conf.householdclust.mod1$y.ellipsepoints) # 1.053940 1.460196

# cluster ascertainment, q1 = 1, give q2 a value
household.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 1, 
                                                   q2 = 0.9, 
                                                   clusters = householdsize))

outputres(household.mod2)
# dispersion_k meanclust        Rc
# 1     1.111628  1.996719 0.9967188
conf.householdclust.mod2 <- conf_ellipse(inputdata = list(a.est = household.mod2$a_n_b.est[1],
                                                          b.est = household.mod2$a_n_b.est[2],
                                                          hessmat = household.mod2$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.householdclust.mod2$x.ellipsepoints) #  1.938963 2.054475
range(conf.householdclust.mod2$y.ellipsepoints) # 0.942947 1.280308

# double ascertainment
household.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.9, 
                                                   q2 = 0.9, 
                                                   clusters = householdsize))

outputres(household.mod3)
# dispersion_k meanclust     Rc
# 1     1.115988    2.1099 1.1099
conf.householdclust.mod3 <- conf_ellipse(inputdata = list(a.est = household.mod3$a_n_b.est[1],
                                                          b.est = household.mod3$a_n_b.est[2],
                                                          hessmat = household.mod3$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.householdclust.mod3$x.ellipsepoints) #  2.045396 2.174405
range(conf.householdclust.mod3$y.ellipsepoints) # 0.9411054 1.2908699

############# carehome size
carehomesize <- clustersize_list$carehome$clustersize

# individual ascertainment, q2 = 1, give q1 a value
carehome.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.9, 
                                                   q2 = 1, 
                                                   clusters = carehomesize))

outputres(carehome.mod1)
# dispersion_k meanclust      Rc
# 1    0.1180123   8.87686 7.87686
conf.carehomeclust.mod1 <- conf_ellipse(inputdata = list(a.est = carehome.mod1$a_n_b.est[1],
                                                          b.est = carehome.mod1$a_n_b.est[2],
                                                          hessmat = carehome.mod1$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.carehomeclust.mod1$x.ellipsepoints) #  -8.435843 26.189562
range(conf.carehomeclust.mod1$y.ellipsepoints) # 0.03545388 0.20057067

# cluster ascertainment, q1 = 1, give q2 a value
carehome.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 1, 
                                                   q2 = 0.7, 
                                                   clusters = carehomesize))

outputres(carehome.mod2)
# dispersion_k meanclust       Rc
#    0.0942729  7.173796 6.173796
conf.carehomeclust.mod2 <- conf_ellipse(inputdata = list(a.est = carehome.mod2$a_n_b.est[1],
                                                          b.est = carehome.mod2$a_n_b.est[2],
                                                          hessmat = carehome.mod2$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.carehomeclust.mod2$x.ellipsepoints) #  -7.463783 21.811374
range(conf.carehomeclust.mod2$y.ellipsepoints) # 0.02766699 0.16087881

# double ascertainment
carehome.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                  inputdata = list(q1 = 0.9, 
                                                   q2 = 0.7, 
                                                   clusters = carehomesize))

outputres(carehome.mod3)
# dispersion_k meanclust       Rc
#    0.0889106  7.320153 6.320153
conf.carehomeclust.mod3 <- conf_ellipse(inputdata = list(a.est = carehome.mod3$a_n_b.est[1],
                                                          b.est = carehome.mod3$a_n_b.est[2],
                                                          hessmat = carehome.mod3$hessian_mat),
                                         conflevel = 0.95,
                                         npoint = 100)
range(conf.carehomeclust.mod3$x.ellipsepoints) #  -7.280224 21.920529
range(conf.carehomeclust.mod3$y.ellipsepoints) # 0.02506432 0.15275688

######## restaurant size
restaurantsize <- allclustcomb.df %>% filter(category == "Restaurants") 
restaurantsize <- restaurantsize$size

# individual ascertainment, q2 = 1, give q1 a value
restaurant.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                 inputdata = list(q1 = 0.8, 
                                                  q2 = 1, 
                                                  clusters = restaurantsize))

outputres(restaurant.mod1)
# dispersion_k meanclust       Rc
#    0.5841888  2.997467 1.997467
conf.restaurantclust.mod1 <- conf_ellipse(inputdata = list(a.est = restaurant.mod1$a_n_b.est[1],
                                                         b.est = restaurant.mod1$a_n_b.est[2],
                                                         hessmat = restaurant.mod1$hessian_mat),
                                        conflevel = 0.95,
                                        npoint = 100)
range(conf.restaurantclust.mod1$x.ellipsepoints) #  2.558494 3.436440
range(conf.restaurantclust.mod1$y.ellipsepoints) # 0.3931208 0.7752567

# cluster ascertainment, q1 = 1, give q2 a value
restaurant.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                 inputdata = list(q1 = 1, 
                                                  q2 = 0.8, 
                                                  clusters = restaurantsize))

outputres(restaurant.mod2)
# dispersion_k meanclust       Rc
#   0.5208021  2.474677 1.474677
conf.restaurantclust.mod2 <- conf_ellipse(inputdata = list(a.est = restaurant.mod2$a_n_b.est[1],
                                                         b.est = restaurant.mod2$a_n_b.est[2],
                                                         hessmat = restaurant.mod2$hessian_mat),
                                        conflevel = 0.95,
                                        npoint = 100)
range(conf.restaurantclust.mod2$x.ellipsepoints) #  2.139410 2.809944
range(conf.restaurantclust.mod2$y.ellipsepoints) # 0.3574778 0.6841264

# double ascertainment
restaurant.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                 inputdata = list(q1 = 0.8, 
                                                  q2 = 0.8, 
                                                  clusters = restaurantsize))

outputres(restaurant.mod3)
# dispersion_k meanclust       Rc
#     0.5025523  2.779084 1.779084
conf.restaurantclust.mod3 <- conf_ellipse(inputdata = list(a.est = restaurant.mod3$a_n_b.est[1],
                                                         b.est = restaurant.mod3$a_n_b.est[2],
                                                         hessmat = restaurant.mod3$hessian_mat),
                                        conflevel = 0.95,
                                        npoint = 100)
range(conf.restaurantclust.mod3$x.ellipsepoints) #  2.374151 3.184018
range(conf.restaurantclust.mod3$y.ellipsepoints) # 0.3370168 0.6680879

###### office
officesize <- clustersize_list$office$clustersize

# individual ascertainment, q2 = 1, give q1 a value
office.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                   inputdata = list(q1 = 0.9, 
                                                    q2 = 1, 
                                                    clusters = officesize))

outputres(office.mod1)
# dispersion_k meanclust       Rc
#    0.3762249  1.445429 0.4454291
conf.officeclust.mod1 <- conf_ellipse(inputdata = list(a.est = office.mod1$a_n_b.est[1],
                                                           b.est = office.mod1$a_n_b.est[2],
                                                           hessmat = office.mod1$hessian_mat),
                                          conflevel = 0.95,
                                          npoint = 100)
range(conf.officeclust.mod1$x.ellipsepoints) # 1.308440 1.582419
range(conf.officeclust.mod1$y.ellipsepoints) # 0.1416088 0.6108409

# cluster ascertainment, q1 = 1, give q2 a value
office.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                   inputdata = list(q1 = 1, 
                                                    q2 = 0.9, 
                                                    clusters = officesize))

outputres(office.mod2)
# dispersion_k meanclust       Rc
#   0.3432065  1.377791 0.377791
conf.officeclust.mod2 <- conf_ellipse(inputdata = list(a.est = office.mod2$a_n_b.est[1],
                                                           b.est = office.mod2$a_n_b.est[2],
                                                           hessmat = office.mod2$hessian_mat),
                                          conflevel = 0.95,
                                          npoint = 100)
range(conf.officeclust.mod2$x.ellipsepoints) #  1.258129 1.497453
range(conf.officeclust.mod2$y.ellipsepoints) # 0.1364052 0.5500077

# double ascertainment
office.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                   inputdata = list(q1 = 0.9, 
                                                    q2 = 0.9, 
                                                    clusters = officesize))

outputres(office.mod3)
# dispersion_k meanclust       Rc
#    0.3370605  1.411617 0.411617
conf.officeclust.mod3 <- conf_ellipse(inputdata = list(a.est = office.mod3$a_n_b.est[1],
                                                           b.est = office.mod3$a_n_b.est[2],
                                                           hessmat = office.mod3$hessian_mat),
                                          conflevel = 0.95,
                                          npoint = 100)
range(conf.officeclust.mod3$x.ellipsepoints) #  1.283280 1.539954
range(conf.officeclust.mod3$y.ellipsepoints) #  0.1304339 0.5436872

####### manual labour
manual_laboursize <- clustersize_list$manual_labour$clustersize

# individual ascertainment, q2 = 1, give q1 a value
manual_labour.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 0.8, 
                                                q2 = 1, 
                                                clusters = manual_laboursize))

outputres(manual_labour.mod1)
# dispersion_k meanclust       Rc
#     0.1367724  2.373062 1.373062
conf.manual_labourclust.mod1 <- conf_ellipse(inputdata = list(a.est = manual_labour.mod1$a_n_b.est[1],
                                                       b.est = manual_labour.mod1$a_n_b.est[2],
                                                       hessmat = manual_labour.mod1$hessian_mat),
                                      conflevel = 0.95,
                                      npoint = 100)
range(conf.manual_labourclust.mod1$x.ellipsepoints) # 1.804436 2.941688
range(conf.manual_labourclust.mod1$y.ellipsepoints) # 0.08316774 0.19037713

# cluster ascertainment, q1 = 1, give q2 a value
manual_labour.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 1, 
                                                q2 = 0.8, 
                                                clusters = manual_laboursize))

outputres(manual_labour.mod2)
# dispersion_k meanclust       Rc
#   0.1255907  2.044567 1.044567
conf.manual_labourclust.mod2 <- conf_ellipse(inputdata = list(a.est = manual_labour.mod2$a_n_b.est[1],
                                                       b.est = manual_labour.mod2$a_n_b.est[2],
                                                       hessmat = manual_labour.mod2$hessian_mat),
                                      conflevel = 0.95,
                                      npoint = 100)
range(conf.manual_labourclust.mod2$x.ellipsepoints) #  1.596644 2.492489
range(conf.manual_labourclust.mod2$y.ellipsepoints) # 0.07751908 0.17366236

# double ascertainment
manual_labour.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 0.8, 
                                                q2 = 0.8, 
                                                clusters = manual_laboursize))

outputres(manual_labour.mod3)
# dispersion_k meanclust       Rc
#    0.1142459  2.156661 1.156661
conf.manual_labourclust.mod3 <- conf_ellipse(inputdata = list(a.est = manual_labour.mod3$a_n_b.est[1],
                                                       b.est = manual_labour.mod3$a_n_b.est[2],
                                                       hessmat = manual_labour.mod3$hessian_mat),
                                      conflevel = 0.95,
                                      npoint = 100)
range(conf.manual_labourclust.mod3$x.ellipsepoints) #  1.671761 2.641562
range(conf.manual_labourclust.mod3$y.ellipsepoints) #  0.06904157 0.15945028

###### nosocomial
nosocomialsize <- clustersize_list$nosocomial$clustersize
# individual ascertainment, q2 = 1, give q1 a value
nosocomial.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                      inputdata = list(q1 = 0.9, 
                                                       q2 = 1, 
                                                       clusters = nosocomialsize))

outputres(nosocomial.mod1)
# dispersion_k meanclust       Rc
#     0.3065185  2.464685 1.464685
conf.nosocomialclust.mod1 <- conf_ellipse(inputdata = list(a.est = nosocomial.mod1$a_n_b.est[1],
                                                              b.est = nosocomial.mod1$a_n_b.est[2],
                                                              hessmat = nosocomial.mod1$hessian_mat),
                                             conflevel = 0.95,
                                             npoint = 100)
range(conf.nosocomialclust.mod1$x.ellipsepoints) # 1.624464 3.304906
range(conf.nosocomialclust.mod1$y.ellipsepoints) # 0.09981972 0.51321720

# cluster ascertainment, q1 = 1, give q2 a value
nosocomial.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                      inputdata = list(q1 = 1, 
                                                       q2 = 0.9, 
                                                       clusters = nosocomialsize))

outputres(nosocomial.mod2)
# dispersion_k meanclust       Rc
#   0.289266  2.278387 1.278387
conf.nosocomialclust.mod2 <- conf_ellipse(inputdata = list(a.est = nosocomial.mod2$a_n_b.est[1],
                                                              b.est = nosocomial.mod2$a_n_b.est[2],
                                                              hessmat = nosocomial.mod2$hessian_mat),
                                             conflevel = 0.95,
                                             npoint = 100)
range(conf.nosocomialclust.mod2$x.ellipsepoints) #  1.524493 3.032280
range(conf.nosocomialclust.mod2$y.ellipsepoints) # 0.09759157 0.48094044

# double ascertainment
nosocomial.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                      inputdata = list(q1 = 0.9, 
                                                       q2 = 0.9, 
                                                       clusters = nosocomialsize))

outputres(nosocomial.mod3)
# dispersion_k meanclust       Rc
#     0.2806839  2.373118 1.373118
conf.nosocomialclust.mod3 <- conf_ellipse(inputdata = list(a.est = nosocomial.mod3$a_n_b.est[1],
                                                              b.est = nosocomial.mod3$a_n_b.est[2],
                                                              hessmat = nosocomial.mod3$hessian_mat),
                                             conflevel = 0.95,
                                             npoint = 100)
range(conf.nosocomialclust.mod3$x.ellipsepoints) #  1.571407 3.174830
range(conf.nosocomialclust.mod3$y.ellipsepoints) #  0.09082425 0.47054346

#### retail and leisure setting
retailsize <- clustersize_list$retail$clustersize
# individual ascertainment, q2 = 1, give q1 a value
retail.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                                   inputdata = list(q1 = 0.7, 
                                                    q2 = 1, 
                                                    clusters = retailsize))

outputres(retail.mod1)
# dispersion_k meanclust       Rc
#    0.05380383  1.873136 0.873136
conf.retailclust.mod1 <- conf_ellipse(inputdata = list(a.est = retail.mod1$a_n_b.est[1],
                                                           b.est = retail.mod1$a_n_b.est[2],
                                                           hessmat = retail.mod1$hessian_mat),
                                          conflevel = 0.95,
                                          npoint = 100)
range(conf.retailclust.mod1$x.ellipsepoints) #  1.019702 2.726569
range(conf.retailclust.mod1$y.ellipsepoints) # 0.01059249 0.09701517

# cluster ascertainment, q1 = 1, give q2 a value
retail.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                                   inputdata = list(q1 = 1, 
                                                    q2 = 0.7, 
                                                    clusters = retailsize))

outputres(retail.mod2)
# dispersion_k meanclust       Rc
#   0.04869164  1.576208 0.5762076
conf.retailclust.mod2 <- conf_ellipse(inputdata = list(a.est = retail.mod2$a_n_b.est[1],
                                                           b.est = retail.mod2$a_n_b.est[2],
                                                           hessmat = retail.mod2$hessian_mat),
                                          conflevel = 0.95,
                                          npoint = 100)
range(conf.retailclust.mod2$x.ellipsepoints) #  0.9840502 2.1683650
range(conf.retailclust.mod2$y.ellipsepoints) # 0.01061871 0.08676456

# double ascertainment
retail.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                                   inputdata = list(q1 = 0.7, 
                                                    q2 = 0.7, 
                                                    clusters = retailsize))

outputres(retail.mod3)
# dispersion_k meanclust       Rc
#    0.04033252  1.635237 0.6352375
conf.retailclust.mod3 <- conf_ellipse(inputdata = list(a.est = retail.mod3$a_n_b.est[1],
                                                           b.est = retail.mod3$a_n_b.est[2],
                                                           hessmat = retail.mod3$hessian_mat),
                                          conflevel = 0.95,
                                          npoint = 100)
range(conf.retailclust.mod3$x.ellipsepoints) #  1.024710 2.245765
range(conf.retailclust.mod3$y.ellipsepoints) #  0.007603704 0.073061331

######## close social indoor
socialsize <- clustersize_list$social$clustersize
# individual ascertainment, q2 = 1, give q1 a value
social.mod1 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 0.7, 
                                                q2 = 1, 
                                                clusters = socialsize))

outputres(social.mod1)
# dispersion_k meanclust       Rc
#    0.1338612   13.9852 12.9852
conf.socialclust.mod1 <- conf_ellipse(inputdata = list(a.est = social.mod1$a_n_b.est[1],
                                                       b.est = social.mod1$a_n_b.est[2],
                                                       hessmat = social.mod1$hessian_mat),
                                      conflevel = 0.95,
                                      npoint = 100)
range(conf.socialclust.mod1$x.ellipsepoints) #  3.473741 24.496668
range(conf.socialclust.mod1$y.ellipsepoints) # 0.06359521 0.20412713

# cluster ascertainment, q1 = 1, give q2 a value
social.mod2 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 1, 
                                                q2 = 0.5, 
                                                clusters = socialsize))

outputres(social.mod2)
# dispersion_k meanclust       Rc
#    0.1025561  8.089177 7.089177
conf.socialclust.mod2 <- conf_ellipse(inputdata = list(a.est = social.mod2$a_n_b.est[1],
                                                       b.est = social.mod2$a_n_b.est[2],
                                                       hessmat = social.mod2$hessian_mat),
                                      conflevel = 0.95,
                                      npoint = 100)
range(conf.socialclust.mod2$x.ellipsepoints) #  2.180466 13.997888
range(conf.socialclust.mod2$y.ellipsepoints) # 0.0478317 0.1572804

# double ascertainment
social.mod3 <- optim_fit_model(init = c(0.5, 0.5), 
                               inputdata = list(q1 = 0.7, 
                                                q2 = 0.5, 
                                                clusters = socialsize))

outputres(social.mod3)
# dispersion_k meanclust       Rc
#    0.0851684  9.170591 8.170591
conf.socialclust.mod3 <- conf_ellipse(inputdata = list(a.est = social.mod3$a_n_b.est[1],
                                                       b.est = social.mod3$a_n_b.est[2],
                                                       hessmat = social.mod3$hessian_mat),
                                      conflevel = 0.95,
                                      npoint = 100)
range(conf.socialclust.mod3$x.ellipsepoints) #  2.389397 15.951784
range(conf.socialclust.mod3$y.ellipsepoints) #  0.03702271 0.13331408

############ visualization

### ellispse
conf.df <- data.frame(
  xconfpoints = c(conf.allclust.mod2$x.ellipsepoints,
                  conf.householdclust.mod2$x.ellipsepoints,
                  conf.restaurantclust.mod2$x.ellipsepoints,
                  conf.socialclust.mod2$x.ellipsepoints,
                  conf.retailclust.mod2$x.ellipsepoints,
                  conf.manual_labourclust.mod2$x.ellipsepoints,
                  conf.officeclust.mod2$x.ellipsepoints,
                  conf.carehomeclust.mod2$x.ellipsepoints,
                  conf.nosocomialclust.mod2$x.ellipsepoints
  ),
  yconfpoints = c(conf.allclust.mod2$y.ellipsepoints,
                  conf.householdclust.mod2$y.ellipsepoints,
                  conf.restaurantclust.mod2$y.ellipsepoints,
                  conf.socialclust.mod2$y.ellipsepoints,
                  conf.retailclust.mod2$y.ellipsepoints,
                  conf.manual_labourclust.mod2$y.ellipsepoints,
                  conf.officeclust.mod2$y.ellipsepoints,
                  conf.carehomeclust.mod2$y.ellipsepoints,
                  conf.nosocomialclust.mod2$y.ellipsepoints
  ),
  cluster = c(rep("All cluster", length(conf.allclust.mod2$x.ellipsepoints)),
              rep("Households", length(conf.householdclust.mod2$x.ellipsepoints)),
              rep("Restaurants", length(conf.restaurantclust.mod2$x.ellipsepoints)),
              rep("Close-social indoor", length(conf.socialclust.mod2$x.ellipsepoints)),
              rep("Retail & leisure", length(conf.retailclust.mod2$x.ellipsepoints)),
              rep("Manual labour", length(conf.manual_labourclust.mod2$x.ellipsepoints)),
              rep("Office work", length(conf.officeclust.mod2$x.ellipsepoints)),
              rep("Care homes", length(conf.carehomeclust.mod2$x.ellipsepoints)),
              rep("Nosocomial", length(conf.nosocomialclust.mod2$x.ellipsepoints)))
  
)
unique(conf.df$cluster)

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
  x = c(outputres(allclust.mod2)[2]%>% as.numeric(), 
        outputres(household.mod2)[2]%>% as.numeric(), 
        outputres(restaurant.mod2)[2]%>% as.numeric(),
        outputres(social.mod2)[2]%>% as.numeric(),
        outputres(retail.mod2)[2]%>% as.numeric(),
        outputres(manual_labour.mod2)[2]%>% as.numeric(),
        outputres(office.mod2)[2]%>% as.numeric(),
        outputres(carehome.mod2)[2]%>% as.numeric(), 
        outputres(nosocomial.mod2)[2]%>% as.numeric()),
  y = c(outputres(allclust.mod2)[1]%>% as.numeric(), 
        outputres(household.mod2)[1]%>% as.numeric(),
        outputres(restaurant.mod2)[1]%>% as.numeric(),
        outputres(social.mod2)[1]%>% as.numeric(),
        outputres(retail.mod2)[1]%>% as.numeric(),
        outputres(manual_labour.mod2)[1]%>% as.numeric(),
        outputres(office.mod2)[1]%>% as.numeric(),
        outputres(carehome.mod2)[1]%>% as.numeric(), 
        outputres(nosocomial.mod2)[1]%>% as.numeric()),
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
pointdf$setting <- pointdf$cluster
conf.df$setting <- conf.df$cluster
  
ellipse.1 <- ggplot() + 
  geom_path(data = conf.df , 
            aes(x = xconfpoints, y = yconfpoints, group = setting, color = setting)) +
  geom_point(data = pointdf , 
             aes(x = x, y = y, group = cluster, color = setting)) + 
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
  scale_x_log10(limits = c(1, 26), breaks = c(1, 1.5, 2, 2.5, 3, 7)) +
#  scale_y_continuous(limits = c(0, 1.5)) +
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
propCI_allclust <- frac_est_CI.2080(confinput = conf.allclust.mod2,
                                    estinput = outputres(allclust.mod2))
propCI_allclust
# $est_fraction
# [1] 0.2592072
# 
# $CI_fraction
# [1] 0.2482493 0.2695629


propCI_household <- frac_est_CI.2080(confinput = conf.householdclust.mod2,
                               estinput = outputres(household.mod2))
propCI_household
# $est_fraction
# [1] 0.3097633
# 
# $CI_fraction
# [1] 0.2938210 0.3228981


propCI_restaurant <- frac_est_CI.2080(confinput = conf.restaurantclust.mod2,
                                estinput = outputres(restaurant.mod2))
propCI_restaurant
# $est_fraction
# [1] 0.2601345

# $CI_fraction
# [1] 0.2142684 0.2952625

propCI_social <- frac_est_CI.2080(confinput = conf.socialclust.mod2,
                                estinput = outputres(social.mod2))
propCI_social
# $est_fraction
# [1] 0.09822056
# 
# $CI_fraction
# [1] 0.04990627 0.13947579

propCI_manual <- frac_est_CI.2080(confinput = conf.manual_labourclust.mod2,
                                estinput = outputres(manual_labour.mod2))
propCI_manual
# $est_fraction
# [1] 0.1041502
# 
# $CI_fraction
# [1] 0.07128696 0.13170386

propCI_office <- frac_est_CI.2080(confinput = conf.officeclust.mod2,
                                 estinput = outputres(office.mod2))
propCI_office
# $est_fraction
# [1] 0.1493419
# 
# $CI_fraction
# [1] 0.09030713 0.18293359

propCI_noso <- frac_est_CI.2080(confinput = conf.nosocomialclust.mod2,
                                estinput = outputres(nosocomial.mod2))
propCI_noso
# $est_fraction
# [1] 0.1864948
# 
# $CI_fraction
# [1] 0.08769571 0.24806173



### for carehome and retail clusters, remove the points where CZ is smaller than 1
ind.exc <- which(conf.carehomeclust.mod2$x.ellipsepoints <= 1)
### note the function is generalized to CZ, just in previous version we termed as R0
### but CZ is equal or larger than R0

propresponsible(R0 = as.numeric(outputres(carehome.mod2)[3]),
                k = as.numeric(outputres(carehome.mod2)[1]),
                prop = 0.8) # 0.09120263

201 - length(ind.exc) # 128
frac_careclust <- numeric(128)
for(i in 1:128){
  xpoint <- conf.carehomeclust.mod2$x.ellipsepoints[-ind.exc][i]
  ypoint <- conf.carehomeclust.mod2$y.ellipsepoints[-ind.exc][i]
  frac_careclust[i] <- propresponsible(R0 = xpoint - 1,
                                       k = ypoint, prop = 0.8)
}
range(frac_careclust)
# 0.0299647 0.1368012

##### do the same for retail cluster
ind.exc <- which(conf.retailclust.mod2$x.ellipsepoints <= 1)
### note the function is generalized to CZ, just in previous version we termed as R0
### but CZ is equal or larger than R0

propresponsible(R0 = as.numeric(outputres(retail.mod2)[3]),
                k = as.numeric(outputres(retail.mod2)[1]),
                prop = 0.8) # 0.04570222

201 - length(ind.exc) # 186
frac_retailclust <- numeric(186)
for(i in 1:186){
  xpoint <- conf.retailclust.mod2$x.ellipsepoints[-ind.exc][i]
  ypoint <- conf.retailclust.mod2$y.ellipsepoints[-ind.exc][i]
  frac_retailclust[i] <- propresponsible(R0 = xpoint - 1,
                                       k = ypoint, prop = 0.8)
}
range(frac_retailclust)
# 0.002012172 0.073982951

prop_20_80_CI_list <- list(
  allclust = propCI_allclust,
  householdclust = propCI_household,
  careclust = range(frac_careclust),
  restaurantclust = propCI_restaurant,
  nosoclust = propCI_noso,
  manualclust = propCI_manual,
  officeclust = propCI_office,
  socialclust = propCI_social,
  retailclust = range(frac_retailclust)
)

# save(prop_20_80_CI_list, file = "prop_20_80_CI_list.rda")


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
  annotate("text", x = pointdf$x[1] + 0.65, y = pointdf$y[1] + 0.03, size = 3.5, label = "All clusters") +
  
  geom_point(data = data.frame(x = pointdf$x[2], y = pointdf$y[2]), size = 1.5) +
  annotate("text", x = pointdf$x[2] + 0.7, y = pointdf$y[2], size = 3.5, label = "Households") +
  
  geom_point(data = data.frame(x = pointdf$x[3], y = pointdf$y[3]), size = 1.5) +
  annotate("text", x = pointdf$x[3] + 0.75, y = pointdf$y[3] - 0.015, size = 3.5, label = "Restaurants") +
  
  geom_point(data = data.frame(x = pointdf$x[4], y = pointdf$y[4]), size = 1.5) +
  annotate("text", x = pointdf$x[4] - 0.5 , y = pointdf$y[4] - 0.02, size = 3.5, label = "Close-social indoor") +
  
  geom_point(data = data.frame(x = pointdf$x[5], y = pointdf$y[5]), size = 1.5) +
  annotate("text", x = pointdf$x[5] + 0.35, y = pointdf$y[5] - 0.01, size = 3.5, label = "Retail & leisure") +
  
  geom_point(data = data.frame(x = pointdf$x[6], y = pointdf$y[6]), size = 1.5) +
  annotate("text", x = pointdf$x[6] + 0.75, y = pointdf$y[6], size = 3.5, label = "Manual labour") +
  
  geom_point(data = data.frame(x = pointdf$x[7], y = pointdf$y[7]), size = 1.5) +
  annotate("text", x = pointdf$x[7] + 0.25, y = pointdf$y[7] - 0.04, size = 3.5, label = "Office work") +
  
  geom_point(data = data.frame(x = pointdf$x[8], y = pointdf$y[8]), size = 1.5) +
  annotate("text", x = pointdf$x[8] + 1, y = pointdf$y[8] + 0.025, size = 3.5, label = "Care homes") +
  
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


pfig3 <- ggarrange(ellipse.1, pcontour, nrow = 1, ncol = 2, align = "h", labels = c("a", "b"))
pfig3


ggsave("Fig3_final.pdf", pfig3, width = 12.5, height = 6.5, units = "in", dpi = 600)




































