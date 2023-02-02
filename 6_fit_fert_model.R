# quantify variation in age-specific fertility based on apoe gene

library(rethinking)
library(rstan)
library(parallel)

# data consists of:
# a matrix of n individuals by y ages, giving a birth event or not, given a year of life (0/1/-99)
# we have an n vector of allele presence or absence (1/0) 
# which we correlate to the prob. of conception

load("d.robj")

# we fit the same model
# first to the binarized 4/ no 4
# then to each level
# provide "n_allele" as data, i.e. whether binary or 3-level score

d$allele_2 <- d$allele
d$allele_3 <- d$allele

# code the bin allele
d$allele_2[d$allele == 33] <- 1
d$allele_2[d$allele == 44|d$allele == 34] <- 2

# code the 3_level allele
d$allele_3[d$allele == 33] <- 1
d$allele_3[d$allele == 34] <- 2
d$allele_3[d$allele == 44] <- 3

# make the NA's an integer for stan
d$con[is.na(d$con)] <- -99

d$bmi[is.na(d$bmi)] <- -99

data <- list()

data[[1]] <- list(n = length(unique(d$id)), 
                  y = ncol(d$con), 
                  allele = as.integer(d$allele_2), 
                  c = d$con, 
                  bmi = scale(d$bmi)[, 1],
                  n_allele = 2)

data[[2]] <- list(n = length(unique(d$id)), 
                  y = ncol(d$con), 
                  allele = as.integer(d$allele_3), 
                  c = d$con, 
                  bmi = scale(d$bmi)[, 1],
                  n_allele = 3) 

m <- stan_model("stan/model_fert.stan")

fit <- mclapply(data, function(data) sampling(m, data = data, chains = 4, cores = 4), 
                mc.cores = 16)

saveRDS(fit, "stanfits/fit-fert.rds")