# quantify variation in proximate determinants of fertility by apoe gene

library(rethinking)
library(rstan)
library(parallel)

# data consists of:
# individual AFR/ALR/IBI/LOSS
# which we correlate to the apoe genotype

load("d.robj")

d$allele_2 <- d$allele

# code the bin allele
d$allele_2[d$allele == 33] <- 1
d$allele_2[d$allele == 44|d$allele == 34] <- 2

# fit the model to ibi

# turn the NAs to somethning for stan
ibi <- d$ibi_b
ibi[is.na(ibi)] <- -99

afr <- d$afr
afr[is.na(afr)] <- -99

alr <- d$alr
alr[is.na(alr)] <- -99

loss <- d$loss_b
loss[is.na(loss)] <- -99
loss <- ifelse(loss == 0, 0, 1)

bmi <- scale(d$bmi)[, 1]
bmi[is.na(bmi)] <- -99

data <- list()

data[[1]] <- list(n = length(d$id), 
                  outcome = ibi, 
                  allele = d$allele_2, 
                  bmi = bmi, 
                  n_allele = 2)

data[[2]] <- data[[3]] <- data[[1]]
data[[2]]$outcome <- afr
data[[3]]$outcome <- alr

m <- stan_model("stan/model_prox.stan")

fit <- lapply(data, function(data) sampling(m, data = data, chains = 4, cores = 4))

saveRDS(fit, "stanfits/fit-prox-2-level.rds")

# fit to the 3 level

# code the 3_level allele
d$allele_3[d$allele == 33] <- 1
d$allele_3[d$allele == 34] <- 2
d$allele_3[d$allele == 44] <- 3

# make data

data[[1]] <- list(n = length(d$id), 
                  outcome = ibi, 
                  allele = d$allele_3, 
                  bmi = bmi, 
                  n_allele = 3)

data[[2]] <- data[[3]] <- data[[1]]
data[[2]]$outcome <- afr
data[[3]]$outcome <- alr

fit <- lapply(data, function(data) sampling(m, data = data, chains = 4, cores = 4))

saveRDS(fit, "stanfits/fit-prox-3-level.rds")

# fit loss model

data <- list(n = length(d$id), 
             loss = loss, 
             allele = d$allele_2, 
             bmi = bmi)

m <- stan_model("stan/model_loss.stan")

fit <- sampling(m, data = data, chains = 4, cores = 4)

saveRDS(fit, "stanfits/fit-loss.rds")
