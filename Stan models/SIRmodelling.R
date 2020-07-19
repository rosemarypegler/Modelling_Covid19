library(rstan)
library(gridExtra)
library(devtools)
library(bayesplot)
library(ggplot2)
library(shinystan)

rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())
set.seed(3)

# imports data and splits the time period
# beta is a time-varying function and needs to be estimated over separate time periods
data <- read.table("italycovid19.csv", header=TRUE, sep=",")
data <- data[1:39,]

date <- data$date
infected <- data$suminfected
removed <- data$sumremoved

# total population
N <- 60e6;

# time period
n_days <- length(infected) 
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions
# these will need to be changed when the starting day is not 0
i0 <- 2
r0 <- 0
s0 <- N - i0  - r0
y0 = c(S = s0, I=i0, R = r0)

# data for Stan
modeldatasir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, infected = infected,
                 removed=removed)

model <- stan_model("finalsir.stan", verbose=TRUE)

# start sampling process 
fitsir <- sampling(model, data = modeldatasir, iter = 5000, warmup = 2000, chains = 4,
                           verbose=TRUE)

# print results from the Stan model
pars=c('beta','R0')
print(fitsir, pars = pars)

# Check whether there were any divergences
rstan::check_divergences(fitsir)

# Check whether the chains are in agreement with each other
stan_dens(fitsir, pars = c("beta"), separate_chains = TRUE)

# plot mean, median, CI and observed values for infected and removed
predictedinfected <- cbind(as.data.frame(summary(
  fitsir, pars = "pred_infected", probs = c(0.05, 0.5, 0.95))$summary), t, infected)
colnames(predictedinfected) <- make.names(colnames(predictedinfected)) # to remove % in the col names

ggplot(predictedinfected, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "blue", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.),size=1) +
  geom_line(mapping = aes(x = t, y = mean)) + 
  geom_point(mapping = aes(y = infected)) +
  labs(x = "Day", y = "Infected")

predictedremoved <- cbind(as.data.frame(summary(
  fitsir, pars = "pred_removed", probs = c(0.05, 0.5, 0.95))$summary), t, removed)
colnames(predictedremoved) <- make.names(colnames(predictedremoved)) # to remove % in the col names

ggplot(predictedremoved, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "pink", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.),size=1) + 
  geom_line(mapping = aes(x = t, y = mean),size=1) +
  geom_point(mapping = aes(y = removed)) +
  labs(x = "Day", y = "Removed")

# Traceplots show whether it has investigated the space effectively
traceplot(fitsir, pars=c("beta"))

# Write results to separate csv files
# This has been done simply to use within Python
df <- data.frame(predinfected)
write.csv(df, 'SIRlockdowninfected.csv')

df <- data.frame(predremoved)
write.csv(df, 'SIRlockdownremoved.csv')

# Another method of looking at the Stan results
# Provides a great deal of information
launch_shinystan(fitsir)
