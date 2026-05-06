library(ergm)

# FITS
#1
m0 <- ergm(g_net ~ edges)
summary(m0)

#2
m1 <- ergm.tapered(g_net ~ edges +
             gwdegree(1.25, fixed = TRUE),
           control = control.ergm.tapered(
             MCMLE.maxit = 60,
             MCMC.samplesize = 15000,
             MCMC.interval = 1000,
             parallel = 4,
             parallel.type = "PSOCK"
           ),
           verbose = TRUE)
summary(m1)
saveRDS(m1, file = "ergm_final_2")

m1 <- readRDS("ergm_final_2")



#3
m3 <- ergm.tapered(g_net ~ edges +
                     gwdegree(1.25, fixed = TRUE) +
                     gwesp(0.1, fixed = TRUE),
                    control = control.ergm.tapered(
                      MCMLE.maxit = 60,
                      MCMC.samplesize = 15000,
                      MCMC.interval = 1000,
                      parallel = 4,
                      parallel.type = "PSOCK"
                    ),
                    verbose = TRUE)
summary(m3)
saveRDS(m3, file = "ergm_final_3.rds")

m3 <- readRDS("ergm_final_3.rds")
#4 

m2 <- ergm.tapered(g_net ~ edges +
                     gwdegree(1.25, fixed = TRUE) +
                     gwesp(0.1, fixed = TRUE) +
                     degree(1:4),
                   control = control.ergm.tapered(
                     MCMLE.maxit = 60,
                     MCMC.samplesize = 15000,
                     MCMC.interval = 1000,
                     parallel = 4,
                     parallel.type = "PSOCK"
                   ),
                   verbose = TRUE)
summary(m2)
saveRDS(m2, file = "ergm_final_4.rds")

m2 <- readRDS("ergm_final_4.rds")

#5 

m4 <- ergm.tapered(g_net ~ edges +
                       gwdegree(1.25, fixed = TRUE) +
                       degree(1:4) +
                      gwnsp(0.5, fixed = TRUE) +
                     nodematch('nucleus') +
                     nodematch('membrane') +
                     nodematch('endomembrane_secretory') +
                     nodematch('complexes') +
                     nodematch('cytoplasm') +
                     nodematch('mitochondrion') +
                     nodematch('other'),
                     control = control.ergm.tapered(
                       MCMLE.maxit = 60,
                       MCMC.samplesize = 25000,
                       MCMC.interval = 1000,
                       parallel = 4,
                       parallel.type = "PSOCK", 
                       MCMC.burnin = 10000
                     ),
                     verbose = TRUE)
summary(m4)

good_m4 <- gof(m4)
plot(good_m4)
good_m4
mcmc.diagnostics(m4)

# reload for when i close R
m4 <- readRDS("ergm_final_5_better_v1.rds")
good_m4 <- readRDS("ergm_final_5_gof_better_v1.rds")

AIC(m0, m1, m2, m3, m4)
AIC(m0) - AIC(m4)


###### m5 and m6 are models not used within the comparison as they proved to be more useful only for prediction

m5 <- ergm.tapered(g_net ~ edges +
                     gwdegree(1.25, fixed = TRUE) +
                     degree(1:4) +
                     gwnsp(0.5, fixed = TRUE) +
                     nodematch('nucleus') +
                     nodematch('membrane') +
                     nodematch('endomembrane_secretory') +
                     nodematch('complexes') +
                     nodematch('cytoplasm') +
                     nodematch('mitochondrion') +
                     nodematch('other'),
                   control = control.ergm.tapered(
                     MCMLE.maxit = 60,
                     MCMC.samplesize = 25000,
                     MCMC.interval = 1000,
                     parallel = 4,
                     parallel.type = "PSOCK", 
                     MCMC.burnin = 10000
                   ),
                   verbose = TRUE)

m6 <- ergm.tapered(g_net ~ edges +
                     gwnsp(2, fixed = TRUE) +
                     gwdegree(1.25, fixed = TRUE) +
                     degree(1:4),
                   control = control.ergm.tapered(
                     MCMLE.maxit = 60,
                     MCMC.samplesize = 25000,
                     MCMC.interval = 1000,
                     parallel = 4,
                     parallel.type = "PSOCK", 
                     MCMC.burnin = 10000
                   ),
                   verbose = TRUE)

summary(m6)


###### goodness of fit tests

good_m1 <- gof(m1)
plot(good_m1)
mcmc.diagnostics(m1)

good_m2 <- gof(m2)
plot(good_m2)
mcmc.diagnostics(m2)

good_m3 <- gof(m3)
plot(good_m3)
mcmc.diagnostics(m3)

good_m4 <- gof(m4)
plot(good_m4)
mcmc.diagnostics(m4)

good_m5 <- gof(m5)
plot(good_m5)
mcmc.diagnostics(m5)
