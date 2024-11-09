library(molMIC)
gent_mics <- example_mics$measurement[example_mics$antibiotic == "gentamicin"]
# gent_mics <- gent_mics[gent_mics > mic_r_breakpoint("E. coli", "gentamicin")]

#' Get the right trough of an MIC distribution
#'
#' @param x A vector of MICs
#'
#' @return The right trough of the MIC distribution (an MIC value)
#'
#' @examples
#' mics <- AMR::as.mic(c(0.25, 4, 8, 8, 0.25, 0.5, 0.125, 1, 1, 2))
#' right_trough(mics) == AMR::as.mic(2)
right_trough <- function(x) {
  x <- AMR::as.mic(x)
  diffs <- x |>
    droplevels() |>
    table() |>
    diff()
  # get first negative diff
  names(diffs[max(which(diffs < 0))])[[1]] |>
    AMR::as.mic()
}

mode_1 <- max(gent_mics)
lvls <- table(droplevels(gent_mics))
right_t <- right_trough(gent_mics)

right <- lvls[AMR::as.mic(names(lvls)) >= right_t]

# samples <- rlnorm(sum(right), meanlog = 0.1, sdlog = 0.5)

grid <- expand.grid(
  meanlog = seq(0.1, 10, length.out = 100),
  sdlog = seq(0.1, 10, length.out = 100)
)

# samples <- vector("list", length = nrow(grid))
# for (i in seq_len(nrow(grid))) {
#   samples[[i]] <- sort(rlnorm(sum(right), meanlog = grid$meanlog[i], sdlog = grid$sdlog[i]))
# }

samples <- apply(grid, 1, \(x) {
  sort(
    do.call(
      rlnorm,
      c(list(n = sum(right)), as.list(x))
    )
  )
}, simplify = FALSE)

max_mic <- max(gent_mics)
# belows <- vector("list", length = length(samples))
# for (i in seq_along(samples)) {
#   belows[[i]] <- samples[[i]][samples[[i]] < max_mic]
# }

compare_against_these <- gent_mics[mic_uncensor(gent_mics) <= max_mic & mic_uncensor(gent_mics) >= right_t] |>
  sort()

rmse <- vector("numeric", length = length(samples))
for (i in seq_along(samples)) {
  sim <- samples[[i]]
  sim <- sim[1:length(compare_against_these)]
  # rmse
  rmse[[i]] <- sqrt(mean((sim - compare_against_these)^2))
}

best <- which.min(rmse)
best_params <- grid[best, ]
best_sims <- samples[[best]]
best_sims <- sample(best_sims)

best_sims_tab <- table(force_mic(best_sims))

gent_mics_mod <- gent_mics
gent_mics_mod[AMR::as.mic(mic_uncensor(gent_mics_mod)) > AMR::as.mic(max_mic)] <- sample(names(best_sims_tab),
                                                                                         size = sum(AMR::as.mic(mic_uncensor(gent_mics_mod)) > AMR::as.mic(max_mic)),
                                                                                         replace = TRUE, prob = best_sims_tab)
# plot gent_mics and gent_mics_mod side by side
par(mfrow = c(1, 2))
plot(gent_mics, main = "Original gent_mics", xlab = "MIC")
plot(mic_uncensor_dens_right(gent_mics), main = "Uncensored", xlab = "MIC")
# plot(mic_uncensor_dens(gent_mics, dist = rnorm, mean = seq(0.1, 10, length.out = 100)), main = "Modified gent_mics (lnorm)", xlab = "MIC")

lvls <- gent_mics %>%
  droplevels() %>%
  levels() %>%
  as.mic()

sequences <- lapply(gent_mics, \(x) {
  growths <- rep(0, sum(lvls < as.mic(x)))
  inhibits <- rep(1, sum(lvls >= as.mic(x)))
  stopifnot(length(growths) + length(inhibits) == length(lvls))
  c(growths, inhibits)
  # seq(from = 1, to = sum(lvls < as.mic(x)))
})

# surv <- lapply(gent_mics, \(x) {
#   if (as.mic(x) == max(lvls)) {
#     return(list(event = 0, time = length(lvls)))
#   }
#   list(
#     event = 1,
#     time = sum(lvls <= as.mic(x))
#   )
# }) %>%
#   # to df
#   do.call(rbind, .)

event <- ifelse(gent_mics == max(lvls), 0, 1)
time <- sapply(gent_mics, \(x) sum(lvls <= as.mic(x)))
data <- data.frame(time = time, event = event)

surv <- Surv(time = time, event = event)
fit <- coxph(Surv(time = time, event = event, type = "right") ~ 1, data = data)
base_haz <- basehaz(fit)

sample_survival_time <- function(base_haz) {
  u <- runif(1)  # Draw a random probability
  # Find the smallest time where cumulative survival <= u
  sampled_time <- min(base_haz$time[exp(-base_haz$hazard) <= u])
  return(sampled_time)
}
cens_data <- data[data$event == 0, ]
cens_data$sampled_time <- sapply(1:nrow(cens_data), function(i) sample_survival_time(base_haz))

# Merge back with original data
final_data <- merge(data, cens_data[, c("time", "sampled_time")], by = "time", all.x = TRUE)
final_data

#====
param_model <- survreg(Surv(time = time, event = event, type="right") ~ 1, data = data, dist = "gaussian")

# Extract censored observations
cens_data <- data[data$event == 0, ]

# Function to sample survival times for censored observations
set.seed(42)  # For reproducibility
sample_parametric_survival_time <- function(model) {
  # Get the scale and shape parameters from the Weibull model
  scale <- model$scale
  shape <- 1 / model$coefficients

  # Generate a random survival time for a censored observation
  # Using the Weibull distribution for survival times
  sampled_time <- rweibull(1, shape = shape, scale = exp(model$coefficients))

  return(sampled_time)
}

# Apply function to each censored observation
cens_data$sampled_time <- sapply(1:nrow(cens_data), function(i) sample_parametric_survival_time(param_model))

# Merge back with the original data
# final_data <- base::merge(data, cens_data[, c("time", "sampled_time")], by = "time", all.x = TRUE)
# final_data
cens_data

#====

# Install packages if not already installed
install.packages("survival")
install.packages("survminer")

# Load the libraries
library(survival)
library(survminer)
library(ggplot2)

filtered_gent <- gent_mics[gent_mics >= trough(gent_mics)]

# Example data
mic_data <- data.frame(
  MIC_value = as.numeric(filtered_gent),  # Measured or censored MIC values
  censoring_status = as.integer(filtered_gent == max(filtered_gent))  # Indicator for censoring status
)

# Creating a survival object
surv_object <- Surv(time = mic_data$MIC_value, event = 1 - mic_data$censoring_status)
# Fit a Kaplan-Meier survival model
km_fit <- survfit(surv_object ~ 1, data = mic_data)
summary(km_fit)
plot(km_fit)
# Fit a Weibull model
weibull_fit <- survreg(surv_object ~ 1, data = mic_data, dist = "weibull")
summary(weibull_fit)
# Predict MIC for censored data using the Weibull model
predicted_MIC <- predict(weibull_fit, type = "quantile", p = 0.5)  # median MIC
mic_data$predicted_MIC <- ifelse(mic_data$censoring_status == 1, predicted_MIC, mic_data$MIC_value)
# Visualize the Kaplan-Meier curve
ggsurvplot(km_fit, data = mic_data, conf.int = TRUE, risk.table = TRUE)

# Draw random samples from the Weibull distribution for censored values
# Use scale and shape parameters from the Weibull fit
shape <- 1 / weibull_fit$scale
scale <- exp(weibull_fit$coefficients)

mic_data$predicted_MIC <- ifelse(
  mic_data$censoring_status == 1,
  rweibull(sum(mic_data$censoring_status), shape = shape, scale = scale),
  mic_data$MIC_value
)

gent_ecoffs <- ecoffs[ecoffs$antibiotic == "AMX",] %>%
  dplyr::select(`0.002`:`512`)
gent_ecoffs_freq <- gent_ecoffs %>%
  t() %>%
  as.numeric()
gent_ecoffs_mics <- colnames(gent_ecoffs)
# sample ecoffs
sampled_ecoffs <- sample(gent_ecoffs_mics, size = sum(gent_ecoffs_freq), replace = TRUE, prob = gent_ecoffs_freq)
sampled_ecoffs <- sample(sampled_ecoffs, size = 1000) %>%
  AMR::as.mic()
censored_ecoffs <- dplyr::if_else(sampled_ecoffs > AMR::as.mic("16"), AMR::as.mic(">16"), sampled_ecoffs)


# Plot the original and modified MIC distributions
# reduce samples for speed
par(mfrow = c(1, 3))
plot(sampled_ecoffs, main = "Original censored_ecoffs", xlab = "MIC")
plot(censored_ecoffs, main = "Modified censored_ecoffs", xlab = "MIC")
plot(mic_uncensor_dens_right(censored_ecoffs), main = "Uncensored", xlab = "MIC")


