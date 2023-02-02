# 1: plot fertility model results

library(rethinking)
library(rstan)
library(parallel)

load("d.robj")

d$allele_2 <- d$allele
d$allele_3 <- d$allele

# code the bin allele
d$allele_2[d$allele == 33] <- 1
d$allele_2[d$allele == 44|d$allele == 34] <- 2

# code the 3_level allele
d$allele_3[d$allele == 33] <- 1
d$allele_3[d$allele == 34] <- 2
d$allele_3[d$allele == 44] <- 3

# load stanfit

fit <- readRDS("stanfits/fit-fert.rds")

### 1: 33 vs. 43/44

# Figure 2: TFR

age <- 1:91

post <- extract.samples(fit[[1]])

p1 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 1]) # APOE3/APOE3
p1 <- inv_logit(p1)
tfr1 <- apply(p1, 1, cumsum)

p2 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 2]) # APOE3/APOE4 or APOE4/APOE4
p2 <- inv_logit(p2)
tfr2 <- apply(p2, 1, cumsum)

png("output/figure_2_tfr.png", 
    res = 250, 
    height = 1400, 
    width = 1800)

plot(NULL, 
     xlim = c(20, 60), 
     ylim = c(1, 12), 
     xlab = "Maternal Age", 
     ylab = "Predicted Number of Children")

lines(apply(tfr1, 1, mean), col = "goldenrod2", lwd = 2)
lines(apply(tfr2, 1, mean), col = "navy", lwd = 2)

shade(apply(tfr1, 1, function(x) HPDI(x, prob = 0.95)), 1:91, col = col.alpha("goldenrod2", 0.2))
shade(apply(tfr2, 1, function(x) HPDI(x, prob = 0.95)), 1:91, col = col.alpha("navy", 0.1))

legend("topleft", 
       col = c("goldenrod2", "navy"), 
       legend = c("APOE3/APOE3", "APOE3/APOE4 | APOE4/APOE4"), 
       lwd = 2)

dev.off()

tfr_contrast <- tfr2[91, ] - tfr1[91, ]
# 0.30 (-0.12, 0.77), p_sign = 0.096 more kids for anyone with e4
stopifnot(abs(mean(tfr_contrast) - 0.3) < 1e-2)
mean(tfr_contrast)
HPDI(tfr_contrast, prob = 0.95)
mean(tfr_contrast < 0)

# Supp Figure 2: prob. birth by age

p1 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 1])
p1 <- inv_logit(p1)
p2 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 2])
p2 <- inv_logit(p2)

png("output/supp_age_specific_prob_birth_2_level.png", 
    res = 250, 
    height = 1400, 
    width = 1800)

plot(NULL, 
     xlim = c(0, 60), 
     ylim = c(0, 0.5), 
     xlab = "age of woman", 
     ylab = "prob. of giving birth")

lines(age, apply(p1, 2, mean), col = "goldenrod2", lwd = 2)
lines(age, apply(p2, 2, mean), col = "navy", lwd = 2)
shade(apply(p1, 2, function(x) HPDI(x, prob = 0.95)), age, col = col.alpha("goldenrod2", 0.2))
shade(apply(p2, 2, function(x) HPDI(x, prob = 0.95)), age, col = col.alpha("navy", 0.1))

o1 <- d$con[d$allele_2 == 1, ]
o1 <- apply(o1, 2, function(x) mean(x, na.rm = TRUE))
points(o1, col = "goldenrod2")

o2 <- d$con[d$allele_2 == 2, ]
o2 <- apply(o2, 2, function(x) mean(x, na.rm = TRUE))
points(o2, col = "navy")

legend("topright", 
       col = c("goldenrod2", "navy"), 
       legend = c("33", "34/44"), 
       lwd = 2, 
       title = "apoe type")

dev.off()

### 2: 33 vs. 43 vs. 44

post <- extract.samples(fit[[2]])

# Figure 2 TFR

p1 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 1]) # 33
p1 <- inv_logit(p1)
tfr1 <- apply(p1, 1, cumsum)

p2 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 2]) # 34
p2 <- inv_logit(p2)
tfr2 <- apply(p2, 1, cumsum)

p3 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 3]) # 44
p3 <- inv_logit(p3)
tfr3 <- apply(p3, 1, cumsum)

tfr_44_33_contrast <- tfr3[91, ] - tfr1[91, ]
# 2.10 (0.65 - 3.63), psign = 0.0035
stopifnot(abs(mean(tfr_44_33_contrast) - 2.18) < 1e-1)
mean(tfr_44_33_contrast)
HPDI(tfr_44_33_contrast, prob = 0.95)
mean(tfr_44_33_contrast < 0)

tfr_34_33_contrast <- tfr2[91, ] - tfr1[91, ]
# 0.15 (-0.30, 0.63) psign = 0.2795
stopifnot(abs(mean(tfr_34_33_contrast) - 0.14) < 1e-1)
mean(tfr_34_33_contrast)
HPDI(tfr_34_33_contrast, prob = 0.95)
mean(tfr_34_33_contrast < 0)

# associated: <0.01, 0.05, 0.19, 0.13, 0.02, 0.19, 0.24, 0.01
# not credibly associated: 0.29
# weakly associated: 0.37, 0.14, 0.24, 0.74
# weak evidence for association: 0.21, 0.33, 0.32
# more strongly associated: 0.02, 0.00, 0.01


png("output/figure_2_tfr_3type.png", 
    res = 250, 
    height = 1400, 
    width = 1800)

plot(NULL, 
     xlim = c(20, 60), 
     ylim = c(1, 12), 
     xlab = "Maternal Age", 
     ylab = "Predicted Number of Children")

lines(apply(tfr1, 1, mean), col = "goldenrod2", lwd = 2)
lines(apply(tfr2, 1, mean), col = "navy", lwd = 2)
lines(apply(tfr3, 1, mean), col = "black", lwd = 2)

shade(apply(tfr1, 1, function(x) HPDI(x, prob = 0.95)), 1:91, col = col.alpha("goldenrod2", 0.2))
shade(apply(tfr2, 1, function(x) HPDI(x, prob = 0.95)), 1:91, col = col.alpha("navy", 0.1))
shade(apply(tfr3, 1, function(x) HPDI(x, prob = 0.95)), 1:91, col = col.alpha("darkgrey", 0.1))

legend("topleft", 
       col = c("goldenrod2", "navy", "black"), 
       legend = c("APOE3/APOE3", "APOE3/APOE4", "APOE4/APOE4"), 
       lwd = 2)

dev.off()

# Supp Figure 2: prob. birth by age

age <- 1:91

p1 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 1])
p1 <- inv_logit(p1)
p2 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 2])
p2 <- inv_logit(p2)
p3 <- sapply(age, function(age) post$mu[, age] + post$b_allele[, 3])
p3 <- inv_logit(p3)

png("output/supp_age_specific_prob_birth_3_level.png", 
    res = 250, 
    height = 1400, 
    width = 1800)

plot(NULL, 
     xlim = c(0, 60), 
     ylim = c(0, 0.5), 
     xlab = "age of woman", 
     ylab = "prob. of giving birth")

lines(age, apply(p1, 2, mean), col = "goldenrod2", lwd = 2)
lines(age, apply(p2, 2, mean), col = "navy", lwd = 2)
lines(age, apply(p3, 2, mean), col = "black", lwd = 2)
shade(apply(p1, 2, function(x) HPDI(x, prob = 0.95)), age, col = col.alpha("goldenrod2", 0.2))
shade(apply(p2, 2, function(x) HPDI(x, prob = 0.95)), age, col = col.alpha("navy", 0.1))
shade(apply(p3, 2, function(x) HPDI(x, prob = 0.95)), age, col = col.alpha("grey", 0.2))

o1 <- d$con[d$allele_3 == 1, ]
o1 <- apply(o1, 2, function(x) mean(x, na.rm = TRUE))
points(o1, col = "goldenrod2")

o2 <- d$con[d$allele_3 == 2, ]
o2 <- apply(o2, 2, function(x) mean(x, na.rm = TRUE))
points(o2, col = "navy")

o3 <- d$con[d$allele_3 == 3, ]
o3 <- apply(o3, 2, function(x) mean(x, na.rm = TRUE))
points(o3, col = "black")

legend("topright", 
       col = c("goldenrod2", "navy", "black"), 
       legend = c("33", "34", "44"), 
       lwd = 2, 
       title = "apoe type")

dev.off()
