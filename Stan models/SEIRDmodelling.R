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
recovered <- data$sumrecovered
deaths <- data$sumdeath

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
e0 <- 2
r0 <- 0
d0 <- 0
s0 <- N - i0  - r0 - d0 - e0
y0 = c(S = s0, E = e0, I=i0, R = r0, D = d0)

# data for Stan
modeldataseird <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, infected = infected,
                      recovered=recovered, deaths=deaths)

model <- stan_model("finalseird.stan", verbose=TRUE)

# start sampling process 
fitseird <- sampling(model, data = modeldataseird, iter = 5000, warmup = 2000, chains = 4,
                    verbose=TRUE)

# print results from the Stan model
pars=c('beta','delta','R0')
print(fitseird, pars = pars)

# Check whether there were any divergences
rstan::check_divergences(fitseird)

# Check whether the chains are in agreement with each other
stan_dens(fitseird, pars = c("beta","delta"), separate_chains = TRUE)

# plot mean, median, CI and observed values for infected and removed
predictedinfected <- cbind(as.data.frame(summary(
  fitseird, pars = "pred_cases", probs = c(0.05, 0.5, 0.95))$summary), t, infected)
colnames(predictedinfected) <- make.names(colnames(predictedinfected)) # to remove % in the col names

ggplot(predictedinfected, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "blue", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.),size=1) +
  geom_line(mapping = aes(x = t, y = mean)) + 
  geom_point(mapping = aes(y = infected)) +
  labs(x = "Day", y = "Infected")

predictedrecovered <- cbind(as.data.frame(summary(
  fitseird, pars = "pred_recovered", probs = c(0.05, 0.5, 0.95))$summary), t, recovered)
colnames(predictedrecovered) <- make.names(colnames(predictedrecovered)) 

ggplot(predictedrecovered, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "pink", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.),size=1) + 
  geom_line(mapping = aes(x = t, y = mean),size=1) +
  geom_point(mapping = aes(y = recovered)) +
  labs(x = "Day", y = "Recovered")

predicteddeaths <- cbind(as.data.frame(summary(
  fitseird, pars = "pred_deaths", probs = c(0.05, 0.5, 0.95))$summary), t, deaths)
colnames(predicteddeaths) <- make.names(colnames(predicteddeaths)) 

ggplot(predicteddeaths, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "red", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.),size=1) + 
  geom_line(mapping = aes(x = t, y = mean),size=1) +
  geom_point(mapping = aes(y = deaths)) +
  labs(x = "Day", y = "Deaths")

# Traceplots show whether it has investigated the space effectively
traceplot(fitseird, pars=c("beta","delta"))

# Write results to separate csv files
# This has been done simply to use within Python
df <- data.frame(predictedinfected)
write.csv(df, 'SEIRDlockdowninfected.csv')

df <- data.frame(predictedrecovered)
write.csv(df, 'SEIRDlockdownrecovered.csv')

df <- data.frame(predicteddeaths)
write.csv(df, 'SEIRDlockdowndeaths.csv')

# Another method of looking at the Stan results
# Provides a great deal of information
launch_shinystan(fitseird)
