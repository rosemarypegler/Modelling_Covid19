library(rstan)
library(gridExtra)
library(devtools)
library(bayesplot)
library(ggplot2)
library(shinystan)

rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
set.seed(3)

# Import the data and split the period
# This is done as beta is a time-varying function and needs to be estimated
# for separate time periods
mydata <- read.table("italycovid19.csv", header=TRUE, sep=",")
mydata <- mydata[1:39,]
mydata

date <- mydata$date
infected <- mydata$suminfected
removed <- mydata$sumremoved

# total population of Italy
N <- 6040e4;

# time period
n_days <- length(infected) 
t <- seq(0, n_days, by = 1)
t0 = 0
t <- t[-1]

#initial conditions
# these will change when the time period does not start from day 0
i0 <- 0
r0 <- 0
s0 <- N - i0  - r0
y0 = c(S = s0, I = i0, R = r0)

# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, infected = infected,
                 removed=removed)

# import the model
model <- stan_model("finalsir.stan", verbose=TRUE)

# run the sampler using 4 chains, 5000 iterations with 2000 warm up iterations
fit_sir_negbin <- sampling(model,
                           data = data_sir,
                           iter = 5000,
                           warmup = 2000,
                           chains = 4,
                           verbose=TRUE)

# Display results for chosen parameters
pars=c('beta','R0')
print(fit_sir_negbin, pars = pars)

# Check whether there were any divergences
rstan::check_divergences(fit_sir_negbin)

# Check whether the chains are in agreement with each other
stan_dens(fit_sir_negbin, pars = c("beta"), separate_chains = TRUE)

# Visualise the predicted numbers
smr_pred <- cbind(as.data.frame(summary(
  fit_sir_negbin, pars = "pred_infected", probs = c(0.05, 0.5, 0.95))$summary), t, infected)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names
color_scheme_set("purple")

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "plum", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = mean)) + 
  geom_point(mapping = aes(y = infected)) +
  labs(x = "Day", y = "Infected")

smr_predremoved <- cbind(as.data.frame(summary(
  fit_sir_negbin, pars = "pred_removed", probs = c(0.05, 0.5, 0.95))$summary), t, removed)
colnames(smr_predremoved) <- make.names(colnames(smr_predremoved))

ggplot(smr_predremoved, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "plum", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.),size=1) + 
  geom_line(mapping = aes(x = t, y = mean),size=1) +
  geom_point(mapping = aes(y = removed)) +
  labs(x = "Day", y = "Removed")

# Traceplots show whether it has investigated the space effectively
traceplot(fit_sir_negbin, pars=c("beta"))

# Write results to separate csv files
# This has been done simply to use within Python
df <- data.frame(smr_pred)
write.csv(df, 'SIRlockdowninfected.csv')

df <- data.frame(smr_predremoved)
write.csv(df, 'SIRlockdownremoved.csv')

# Another method of looking at the Stan results
# Provides a great deal of information
launch_shinystan(fit_sir_negbin)
