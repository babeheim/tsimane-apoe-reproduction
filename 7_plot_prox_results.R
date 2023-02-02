# quantify variation in proximate determinants of fertility by apoe gene

library(rethinking)
library(rstan)
library(parallel)

load("d.robj")

# first to the binarized 4/ no 4
# then to each level
# provide "n_allele" as data, i.e. whether binary or 3-level score

d$allele_2 <- d$allele

# code the bin allele
d$allele_2[d$allele == 33] <- 1
d$allele_2[d$allele == 44|d$allele == 34] <- 2

# load fit
fit <- readRDS("stanfits/fit-prox-2-level.rds")

# plot IBI

post <- extract.samples(fit[[1]])

mu1 <- post$a + post$b_allele[, 1] # 33
mu2 <- post$a + post$b_allele[, 2] # 34 or 44

png("output/ibi_cat.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(0, 3), 
     ylim = c(1.5, 3), 
     ylab = "Predicted Interbirth Interval (years)", 
     xlab = "", 
     xaxt = "n")

points(1, mean(mu1), cex = 2, lwd = 2)
points(2, mean(mu2), cex = 2, lwd = 2)
arrows(x0 = 1,
       x1 = 1, 
       y0 = HPDI(mu1, prob = 0.95)[1], 
       y1 = HPDI(mu1, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)
arrows(x0 = 2,
       x1 = 2, 
       y0 = HPDI(mu2, prob = 0.95)[1], 
       y1 = HPDI(mu2, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)

mtext(side = 1, 
      at = c(1, 2),
      c("APOE3/APOE3", "APOE4+"))

dev.off()

png("output/ibi_post_dens.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(1, 3), 
     ylim = c(0, 12), 
     yaxt = "n", 
     ylab = "", 
     xlab = "Predicted Interbirth Interval")

polygon(density(mu1), 
        lwd = 2, 
        border = col.alpha("goldenrod2", 0.4), 
        col = col.alpha("goldenrod2", 0.4))

polygon(density(mu2), 
        lwd = 2, 
        border = col.alpha("navy", 0.4), 
        col = col.alpha("navy", 0.4))

abline(v = mean(mu1), col = col.alpha("goldenrod2", 0.7), lty = 2, lwd = 2)
abline(v = mean(mu2), col = col.alpha("navy", 0.7), lty = 2, lwd = 2)

legend("topleft",
       legend = c("APOE3/APOE3", "APOE4+"),
       fill = c(col.alpha("goldenrod2", 0.5),
                col.alpha("navy", 0.4)), 
       border = c(col.alpha("goldenrod2", 0.5),
                  col.alpha("navy", 0.4)))

dev.off()

# get means 

ibi_contrast <- (mu1 - mu2) # 33 vs 34/44
# 0.24 (0.06 - 0.41), psign = 0.004
stopifnot(abs(mean(ibi_contrast) - (0.24)) < 1e-1)
mean(ibi_contrast)
HPDI(ibi_contrast, prob = 0.95)
mean(ibi_contrast < 0)

ibi_perc_change <- 100 * (mu1 - mu2) / mu1
# 10.3 (2.87-17.95), psign = 0.004
mean(ibi_perc_change)
HPDI(ibi_perc_change, prob = 0.95)
mean(ibi_perc_change < 0) 

# plot AFR

post <- extract.samples(fit[[2]])

mu1 <- post$a + post$b_allele[, 1] # 33
mu2 <- post$a + post$b_allele[, 2] # 34 or 44

png("output/afr_cat.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(0, 3), 
     ylim = c(15, 19), 
     ylab = "Predicted Age at First Reproduction", 
     xlab = "", 
     xaxt = "n")

points(1, mean(mu1), cex = 2, lwd = 2)
points(2, mean(mu2), cex = 2, lwd = 2)
arrows(x0 = 1,
       x1 = 1, 
       y0 = HPDI(mu1, prob = 0.95)[1], 
       y1 = HPDI(mu1, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)
arrows(x0 = 2,
       x1 = 2, 
       y0 = HPDI(mu2, prob = 0.95)[1], 
       y1 = HPDI(mu2, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)

mtext(side = 1, 
      at = c(1, 2),
      c("APOE3/APOE3", "APOE4+"))

dev.off()

png("output/afr_post_dens.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(15, 20), 
     ylim = c(0, 4), 
     yaxt = "n", 
     ylab = "", 
     xlab = "Predicted Age at First Reproduction")

polygon(density(mu1), 
        lwd = 2, 
        border = col.alpha("goldenrod2", 0.4), 
        col = col.alpha("goldenrod2", 0.4))

polygon(density(mu2), 
        lwd = 2, 
        border = col.alpha("navy", 0.4), 
        col = col.alpha("navy", 0.4))

abline(v = mean(mu1), col = col.alpha("goldenrod2", 0.7), lty = 2, lwd = 2)
abline(v = mean(mu2), col = col.alpha("navy", 0.7), lty = 2, lwd = 2)

legend("topleft",
       legend = c("APOE3/APOE3", "APOE4+"),
       fill = c(col.alpha("goldenrod2", 0.5),
                col.alpha("navy", 0.4)), 
       border = c(col.alpha("goldenrod2", 0.5),
                  col.alpha("navy", 0.4)))

dev.off()

# get means

afr_contrast <- mu1 - mu2
# 0.79 (0.17 - 1.45), psign = 0.007
stopifnot(abs(mean(afr_contrast) - 0.8) < 1e-1)
mean(afr_contrast)
HPDI(afr_contrast, prob = 0.95)
mean(afr_contrast < 0)

# plot ALR

post <- extract.samples(fit[[3]])

mu1 <- post$a + post$b_allele[, 1] # 33
mu2 <- post$a + post$b_allele[, 2] # 34 or 44

png("output/alr_cat.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(0, 3), 
     ylim = c(30, 40), 
     ylab = "Predicted Age at Last Reproduction", 
     xlab = "", 
     xaxt = "n")

points(1, mean(mu1), cex = 2, lwd = 2)
points(2, mean(mu2), cex = 2, lwd = 2)
arrows(x0 = 1,
       x1 = 1, 
       y0 = HPDI(mu1, prob = 0.95)[1], 
       y1 = HPDI(mu1, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)
arrows(x0 = 2,
       x1 = 2, 
       y0 = HPDI(mu2, prob = 0.95)[1], 
       y1 = HPDI(mu2, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)

mtext(side = 1, 
      at = c(1, 2),
      c("APOE3/APOE3", "APOE4+"))

dev.off()

png("output/alr_post_dens.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(28, 38), 
     ylim = c(0, 4), 
     yaxt = "n", 
     ylab = "", 
     xlab = "Predicted Age at Last Reproduction")

polygon(density(mu1), 
        lwd = 2, 
        border = col.alpha("goldenrod2", 0.4), 
        col = col.alpha("goldenrod2", 0.4))

polygon(density(mu2), 
        lwd = 2, 
        border = col.alpha("navy", 0.4), 
        col = col.alpha("navy", 0.4))

abline(v = mean(mu1), col = col.alpha("goldenrod2", 0.7), lty = 2, lwd = 2)
abline(v = mean(mu2), col = col.alpha("navy", 0.7), lty = 2, lwd = 2)

legend("topleft",
       legend = c("APOE3/APOE3", "APOE4+"),
       fill = c(col.alpha("goldenrod2", 0.5),
                col.alpha("navy", 0.4)), 
       border = c(col.alpha("goldenrod2", 0.5),
                  col.alpha("navy", 0.4)))

dev.off()

# get means ALR

mean(mu1)
HPDI(mu1, prob = 0.95)
mean(mu2)
HPDI(mu2, prob = 0.95)

alr_contrast <- mu1 - mu2 # 33 minus 34/44
# 3.94 (2.62, 5.37) p_sign < 0.0001
mean(alr_contrast)
HPDI(alr_contrast, prob = 0.95)
mean(alr_contrast < 0)

# plot loss results

fit <- readRDS("stanfits/fit-loss.rds")

post <- extract.samples(fit)

p1 <- inv_logit(post$a + post$b_allele[, 1])
p2 <- inv_logit(post$a + post$b_allele[, 2])

png("output/loss_cat.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(0, 3), 
     ylim = c(0, 1), 
     ylab = "Predicted Foetal Loss Probability", 
     xlab = "", 
     xaxt = "n")

points(1, mean(p1), cex = 2, lwd = 2)
points(2, mean(p2), cex = 2, lwd = 2)
arrows(x0 = 1,
       x1 = 1, 
       y0 = HPDI(p1, prob = 0.95)[1], 
       y1 = HPDI(p1, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)
arrows(x0 = 2,
       x1 = 2, 
       y0 = HPDI(p2, prob = 0.95)[1], 
       y1 = HPDI(p2, prob = 0.95)[2], 
       length = 0, 
       lwd = 2)

mtext(side = 1, 
      at = c(1, 2),
      c("APOE3/APOE3", "APOE4+"))

dev.off()

png("output/loss_post_dens.png", 
    res = 250, 
    height = 1400, 
    width = 1600)

plot(NULL, 
     xlim = c(0, 1), 
     ylim = c(0, 30),
     yaxt = "n", 
     ylab = "", 
     xlab = "Predicted Foetal Loss Probability")

polygon(density(p1), 
        lwd = 2, 
        border = col.alpha("goldenrod2", 0.4), 
        col = col.alpha("goldenrod2", 0.4))

polygon(density(p2), 
        lwd = 2, 
        border = col.alpha("navy", 0.4), 
        col = col.alpha("navy", 0.4))

abline(v = mean(p1), col = col.alpha("goldenrod2", 0.7), lty = 2, lwd = 2)
abline(v = mean(p2), col = col.alpha("navy", 0.7), lty = 2, lwd = 2)

legend("topleft",
       legend = c("APOE3/APOE3", "APOE4+"),
       fill = c(col.alpha("goldenrod2", 0.5),
                col.alpha("navy", 0.4)), 
       border = c(col.alpha("goldenrod2", 0.5),
                  col.alpha("navy", 0.4)))

dev.off()

# calculate odds ratio

p1 <- post$a + post$b_allele[, 1]
p2 <- post$a + post$b_allele[, 2]

odds_ratio <- exp(p2) / exp(p1)
# 0.72 (0.49 - 0.96), psign = 0.0205
abs(mean(odds_ratio) - 0.7191611) < 1e-5
mean(odds_ratio)
HPDI(odds_ratio, prob = 0.95)
mean(odds_ratio > 1) 





# compute all the 3-level results

# load fit
fit <- readRDS("stanfits/fit-prox-3-level.rds")

# IBI

post <- extract.samples(fit[[1]])

mu1 <- post$a + post$b_allele[, 1] # 33
mu2 <- post$a + post$b_allele[, 2] # 34
mu3 <- post$a + post$b_allele[, 3] # 44

# get means 

ibi_33_34_contrast <- mu1 - mu2
# 0.23 (0.05-0.42) p_sign < 0.008
mean(ibi_33_34_contrast)
HPDI(ibi_33_34_contrast, prob = 0.95)
mean(ibi_33_34_contrast < 0)

ibi_33_44_contrast <- mu1 - mu3
# 0.46 (-0.28 - 1.15) p_sign = 0.106
mean(ibi_33_44_contrast)
HPDI(ibi_33_44_contrast, prob = 0.95)
mean(ibi_33_44_contrast < 0)

# AFR

# get means

post <- extract.samples(fit[[2]])

mu1 <- post$a + post$b_allele[, 1]
mu2 <- post$a + post$b_allele[, 2]
mu3 <- post$a + post$b_allele[, 3]

afr_33_34_contrast <- mu1 - mu2
# 0.65 (-0.04, 1.29), p_sign = 0.030
mean(afr_33_34_contrast)
HPDI(afr_33_34_contrast, prob = 0.95)
mean(afr_33_34_contrast < 0)

afr_33_44_contrast <- mu1 - mu3
# 3.96 (2.40, 5.63) p_sign < 0.001
mean(afr_33_44_contrast)
HPDI(afr_33_44_contrast, prob = 0.95)
mean(afr_33_44_contrast < 0)
