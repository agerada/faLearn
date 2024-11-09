censor <- "16"

gent_ecoffs <- ecoffs[ecoffs$antibiotic == "GEN",] %>%
  dplyr::select(`0.002`:`512`)
gent_ecoffs_freq <- gent_ecoffs %>%
  t() %>%
  as.numeric()
gent_ecoffs_mics <- colnames(gent_ecoffs)
# sample ecoffs
sampled_ecoffs <- sample(gent_ecoffs_mics, size = sum(gent_ecoffs_freq), replace = TRUE, prob = gent_ecoffs_freq)
sampled_ecoffs <- sample(sampled_ecoffs, size = 2000) %>%
  AMR::as.mic()
censored_ecoffs <- dplyr::if_else(sampled_ecoffs > AMR::as.mic(censor), AMR::as.mic(paste0(">", censor)), sampled_ecoffs)

# calculate lambda
lambda <- AMR::proportion_SI(AMR::as.sir(censored_ecoffs,
                                        mo = "E. coli",
                                        ab = "GEN"))
n <- 2000
sweep_n <- 10

# Generate a mixture of two gaussians, with mixture lambda
x <- c(rnorm(n * lambda, mean = -2, sd = 1), rnorm(n * (1 - lambda), mean = 3, sd = 1))
hist(x)

bp <-  log2(as.numeric(mic_s_breakpoint("E. coli", "GEN", accept_ecoff = T)))

grid <- expand.grid(
  mean1 = seq(-8, bp, length.out = sweep_n),
  sd1 = seq(0.5, 3, length.out = sweep_n),
  mean2 = seq(bp, 10, length.out = sweep_n),
  sd2 = seq(0.5, 3, length.out = sweep_n)
)

# iterate over grid1 and grid2, generating samples from a gaussian mixture
# with parameters for left distribution from grid1 and right distribution from grid2.
samples <- lapply(1:nrow(grid), \(i) {
  left <- rnorm(n * lambda, mean = grid$mean1[i], sd = grid$sd1[i])
  right <- rnorm(n * (1 - lambda), mean = grid$mean2[i], sd = grid$sd2[i])
  c(left, right)
})

log2_censored <- mic_uncensor(censored_ecoffs, method = "scale") %>% log2

# iterate over samples, measure the rmse compared to log2_censored, but only
# for values <= censor
log2_censored_compare <- log2_censored[log2_censored <= log2(as.numeric(censor))]
rmses <- sapply(samples, \(sample) {
  sample <- sample[sample <= log2(as.numeric(censor))]
  boot <- sample(sample, size = length(log2_censored_compare), replace = TRUE)
  sqrt(mean((boot - log2_censored_compare)^2))
})

best_sample <- samples[[which.min(rmses)]]
par(mfrow = c(1, 3))
hist(log2(sampled_ecoffs), main = "Original", xlab = "log2(MIC)")
hist(log2_censored, main = "Censored log2", xlab = "log2(MIC)")
hist(best_sample, main = "Best sample", xlab = "log2(MIC)")

# Plot the original and modified MIC distributions
# reduce samples for speed
par(mfrow = c(1, 3))
plot(sampled_ecoffs, main = "Original censored_ecoffs", xlab = "MIC")
plot(censored_ecoffs, main = "Modified censored_ecoffs", xlab = "MIC")
plot(mic_u(censored_ecoffs), main = "Uncensored", xlab = "MIC")
