# Project 1 - Weibull model for the lung cancer data

# install.packages(c("rjags", "coda", "survival"))  # only needed once

library(rjags)
library(coda)
library(survival)

dir.create("outputs/tables", recursive = TRUE)
dir.create("outputs/figures", recursive = TRUE)

parameters <- c("beta0", "beta_trt", "beta_karno",
                "beta_small", "beta_adeno", "beta_large", "k")

# Data preparation

dat <- subset(veteran, status == 1)
N   <- nrow(dat)

dat$trt2 <- as.numeric(dat$trt == 2)
dat$celltype <- factor(dat$celltype,
                       levels = c("squamous", "smallcell", "adeno", "large"))
dat$smallcell <- as.numeric(dat$celltype == "smallcell")
dat$adeno     <- as.numeric(dat$celltype == "adeno")
dat$large     <- as.numeric(dat$celltype == "large")

karno_mean  <- mean(dat$karno)
dat$karno_c <- dat$karno - karno_mean

cat("Data: N =", N, ", mean(karno) =", round(karno_mean, 2), "\n")
print(table(dat$celltype))

my.data <- list(
  N = N, time = dat$time, trt2 = dat$trt2, karno_c = dat$karno_c,
  smallcell = dat$smallcell, adeno = dat$adeno, large = dat$large
)

# Q1 is the BUGS model file. Here I run the MCMC and check convergence.

set.seed(2026)
my.inits <- list(
  list(beta0 = 0, beta_trt = 0, beta_karno = 0, beta_small = 0,
       beta_adeno = 0, beta_large = 0, k = 1,
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 101),
  list(beta0 = 1, beta_trt = 0.5, beta_karno = 0.01, beta_small = -0.5,
       beta_adeno = -1, beta_large = 0.2, k = 0.8,
       .RNG.name = "base::Marsaglia-Multicarry", .RNG.seed = 202),
  list(beta0 = -1, beta_trt = -0.5, beta_karno = -0.01, beta_small = 0.5,
       beta_adeno = 1, beta_large = -0.2, k = 1.2,
       .RNG.name = "base::Super-Duper", .RNG.seed = 303)
)

jags_fit <- jags.model("model/project1_weibull.bug",
                       data = my.data, inits = my.inits,
                       n.chains = 3, n.adapt = 1000)
update(jags_fit, 10000)
weibull_sim <- coda.samples(jags_fit, parameters, n.iter = 20000, thin = 1)

samps <- as.matrix(weibull_sim)

# Some basic convergence diagnostics for the report
gr  <- gelman.diag(weibull_sim, multivariate = TRUE, autoburnin = FALSE)
ess <- effectiveSize(weibull_sim)
gw  <- geweke.diag(weibull_sim)
geweke_max <- sapply(parameters, function(p) {
  max(abs(sapply(gw, function(g) g$z[p])))
})
diag_tab <- data.frame(
  parameter = parameters,
  Rhat = round(gr$psrf[parameters, "Point est."], 4),
  ESS = round(ess[parameters]),
  Geweke_maxZ = round(geweke_max, 4)
)
write.csv(diag_tab, "outputs/tables/q2_diagnostics.csv", row.names = FALSE)
print(diag_tab)

pdf("outputs/figures/q2_traceplots.pdf", width = 10, height = 8)
plot(weibull_sim)
dev.off()

# Small helper functions used below. In this parameterisation:
# log(lambda) = -eta and S(t) = exp(-lambda * t^k).

eta_lin <- function(trt2, karno_c, smallcell, adeno, large) {
  samps[, "beta0"] +
    samps[, "beta_trt"]   * trt2 +
    samps[, "beta_karno"] * karno_c +
    samps[, "beta_small"] * smallcell +
    samps[, "beta_adeno"] * adeno +
    samps[, "beta_large"] * large
}

post_sum <- function(x) {
  data.frame(
    mean = round(mean(x), 4),
    sd = round(sd(x), 4),
    median = round(median(x), 4),
    q2.5 = round(quantile(x, 0.025), 4),
    q97.5 = round(quantile(x, 0.975), 4),
    row.names = NULL
  )
}

# Q3: posterior summaries

st <- summary(weibull_sim)$statistics
qu <- summary(weibull_sim)$quantiles
q3_tab <- data.frame(
  parameter = parameters,
  mean = round(st[parameters, "Mean"], 4),
  sd = round(st[parameters, "SD"], 4),
  median = round(qu[parameters, "50%"], 4),
  q2.5 = round(qu[parameters, "2.5%"], 4),
  q97.5 = round(qu[parameters, "97.5%"], 4)
)
write.csv(q3_tab, "outputs/tables/q3_posterior_summary.csv", row.names = FALSE)

# Q4: hazard ratios, using HR = exp(-beta)

hr_names <- c("beta_trt", "beta_karno", "beta_small", "beta_adeno", "beta_large")
hr_tab <- NULL
for (b in hr_names) {
  x <- exp(-samps[, b])
  row <- cbind(covariate = b, post_sum(x))
  hr_tab <- rbind(hr_tab, row)
}
write.csv(hr_tab, "outputs/tables/q4_hazard_ratios.csv", row.names = FALSE)

pdf("outputs/figures/q4_hazard_ratio_forest.pdf", width = 8, height = 5)
par(mar = c(4, 10, 3, 2))
y <- nrow(hr_tab):1
plot(hr_tab$mean, y, xlim = range(c(hr_tab$q2.5, hr_tab$q97.5, 1)),
     pch = 19, yaxt = "n", xlab = "Hazard ratio", ylab = "",
     main = "Hazard ratios (95% credible intervals)")
segments(hr_tab$q2.5, y, hr_tab$q97.5, y)
abline(v = 1, lty = 2)
axis(2, at = y, labels = hr_names, las = 1)
dev.off()

# Q5: median survival by cell type, for treatment 1

cell_types <- c("squamous", "smallcell", "adeno", "large")
karno_by_type <- tapply(dat$karno, dat$celltype, mean)

q5_tab <- NULL
for (ct in cell_types) {
  sc <- as.numeric(ct == "smallcell")
  ad <- as.numeric(ct == "adeno")
  lg <- as.numeric(ct == "large")
  kc <- karno_by_type[[ct]] - karno_mean
  eta <- eta_lin(0, kc, sc, ad, lg)
  med <- (log(2) / exp(-eta))^(1 / samps[, "k"])
  row <- cbind(data.frame(celltype = ct, trt = 1), post_sum(med))
  q5_tab <- rbind(q5_tab, row)
}
write.csv(q5_tab, "outputs/tables/q5_median_survival_by_celltype.csv", row.names = FALSE)

pdf("outputs/figures/q5_median_survival_caterpillar.pdf", width = 8, height = 5)
par(mar = c(4, 10, 3, 2))
y <- nrow(q5_tab):1
plot(q5_tab$mean, y, xlim = range(c(q5_tab$q2.5, q5_tab$q97.5)),
     pch = 19, yaxt = "n", xlab = "Median survival (days)", ylab = "")
segments(q5_tab$q2.5, y, q5_tab$q97.5, y)
axis(2, at = y, labels = q5_tab$celltype, las = 1)
dev.off()

# Q6: survival at day 60 for the two profiles

eta_A <- eta_lin(0, 0, 0, 0, 0)
eta_B <- eta_lin(1, 0, 0, 0, 0)
SA <- exp(-exp(-eta_A) * 60^samps[, "k"])
SB <- exp(-exp(-eta_B) * 60^samps[, "k"])

q6_tab <- rbind(
  cbind(profile = "A", post_sum(SA)),
  cbind(profile = "B", post_sum(SB))
)
write.csv(q6_tab, "outputs/tables/q6_survival_at_60.csv", row.names = FALSE)

pdf("outputs/figures/q6_survival_at_60_density.pdf", width = 8, height = 5)
plot(density(SA), col = "blue", lwd = 2, xlab = "S(60)", ylab = "Density",
     main = "Profile A vs B")
lines(density(SB), col = "red", lwd = 2)
legend("topright", c("A (trt 1)", "B (trt 2)"), col = c("blue", "red"), lwd = 2)
dev.off()

# Q7: posterior probability that profile A has larger S(60)

prob_q7 <- mean(SA > SB)
q7_tab <- data.frame(
  estimate = round(prob_q7, 4),
  n_draws = length(SA)
)
write.csv(q7_tab, "outputs/tables/q7_prob_SA_gt_SB.csv", row.names = FALSE)

# Q8: probability that the median survival is above 100 days

q8_tab <- NULL
for (ct in cell_types) {
  sc <- as.numeric(ct == "smallcell")
  ad <- as.numeric(ct == "adeno")
  lg <- as.numeric(ct == "large")
  kc <- karno_by_type[[ct]] - karno_mean
  eta <- eta_lin(0, kc, sc, ad, lg)
  med <- (log(2) / exp(-eta))^(1 / samps[, "k"])
  row <- data.frame(
    celltype = ct,
    prob_median_gt_100 = round(mean(med > 100), 4)
  )
  q8_tab <- rbind(q8_tab, row)
}
write.csv(q8_tab, "outputs/tables/q8_prob_median_gt_100.csv", row.names = FALSE)

cat("\nDone. Tables in outputs/tables/, figures in outputs/figures/\n")
