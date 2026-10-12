# Loading and inspecting the data ************************************************************************
data <- read.table("transformed_power_data.csv", header=T, sep = ",")


Sys.setlocale("LC_TIME", "English")
data$start_datetime <- as.POSIXct(data$start_datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")  
# tz="UTC" is necessary to avoid artificial duplicates caused by the Spring DST transition

head(data)
tail(data)
print(data[35927:35934,])
unique(data$notes) # "weekday"        "weekend"        "vacation"       "COVID_lockdown"
hist(data$month_of_year) # june a bit overdistributed since the dataset starts with june and ends with early july
unique(data$month_of_year)
hist(data$day_of_week) # just checking, seems fine
hist(data$hour_of_day) # too many zeroes?
unique(data$hour_of_day) # only the correct values
table(data$hour_of_day) # nope its just histogram putting 0 and 1 in the same column

summary(data)

# ************* Data looks OK

# NA check
length(data$value_kWh) == length(na.contiguous(data$value_kWh))  # == 0 [OK]
# ************* no NA found, fine

# Duplicates check
nrow(data) - nrow(unique(data))  # == 0 [OK]

length(data$start_datetime) - length(unique(data$start_datetime))  # == 0 [OK]

# data[duplicated(data$start_datetime), ]

# # show all duplicated rows
# data[duplicated(data$start_datetime) | duplicated(data$start_datetime, fromLast = TRUE), ]  # strange behavior: see hour in "start_datetime" and "hour_of_day"

# ************* Duplicates not found

# **************************************************************************************************
## Plotting the target variable
hist(data$value_kWh, prob=T, breaks = 50)
plot(data$start_datetime, data$value_kWh, type = "l") 

## Inspecting the relationships of power consumption to other variables ***************************
# since these are not continuous, correlation coefficient or scatter plots cannot be used
# boxplots show the relationships nicely

boxplot(value_kWh ~ day_of_week, data = data)
boxplot(value_kWh ~ month_of_year, data = data)
boxplot(value_kWh ~ hour_of_day, data = data)
boxplot(value_kWh ~ notes, data = data)




# **************************************************************************************************
# **************************************************************************************************



### 3-state model - with self-written likelihood function
## We will be using gamma distribution - positive right skewed data

# for a 3-state HMM with state-dependent gamma distributions of the power consumption,
# starting values could be:
mu0 <- c(0.3, 0.9, 1.5)
sigma0 <- c(2, 4, 5)
etas0 <- rep(0, 6) # transition matrix parameters starting values, 6 times 0
theta <- c(etas0, mu0, sigma0)

# transform the starting values for theta.star
theta.star <- c(theta[1:6], log(theta[7:12])) # etas are unconstrained

transform_3x3 <- function(theta.star) {
  # first row
  eta1 <- theta.star[1]
  eta2 <- theta.star[2]
  den <- 1 + exp(eta1) + exp(eta2) # exp so they are positive, denominator to make it all sum to 1
  g11 <- exp(eta1)/den
  g12 <- exp(eta2)/den
  g13 <- 1/den
  
  eta1 <- theta.star[3]
  eta2 <- theta.star[4]
  den <- 1 + exp(eta1) + exp(eta2) 
  g21 <- exp(eta1)/den
  g22 <- exp(eta2)/den
  g23 <- 1/den
  
  eta1 <- theta.star[5]
  eta2 <- theta.star[6]
  den <- 1 + exp(eta1) + exp(eta2) 
  g31 <- exp(eta1)/den
  g32 <- exp(eta2)/den
  g33 <- 1/den
  
  Gamma <- matrix(
    c(g11,g12,g13,
      g21,g22,g23,
      g31,g32,g33),
    nrow=3, byrow=TRUE)
  
  return(Gamma)
}


## ************ Likelihood function dealing with parameter constraints
# and returning the negative log-likelihood
# phi is there to deal with numerical underflow
# allprobs to deal with missing values. (But actually there are none in our dataset, so this is not necessary)

llk <- function(theta.star, x) { # unconstrained parameters theta.star
  Gamma <- transform_3x3(theta.star[1:6]) # retransform probabilities - to be in [0, 1], each rowsum 1
  # retransform mean and standard deviations - to be positive
  mu <- exp(theta.star[7:9])
  sigma <- exp(theta.star[10:12])
  delta <- rep(1 / 3, 3) # 3-states initial probs
  
  allprobs <- matrix(1, length(x), 3) # allprobs = matrix full of 1s, for every observation
  # for t's when allprobs is 1 the diagonal matrix we create in the 
  # forward algorithm becomes the identity matrix
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind( # now only put density values for non-missing observations
    dgamma(x[ind], shape = mu[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu[1]),
    dgamma(x[ind], shape = mu[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu[2]),
    dgamma(x[ind], shape = mu[3]^2 / sigma[3]^2, scale = sigma[3]^2 / mu[3])
  )
  foo <- delta %*% diag(allprobs[1, ]) # alpha_1
  l <- log(sum(foo)) # this will add up to the log-likelihood
  phi <- foo / sum(foo) # standardized alpha_1: phi_1
  for (t in 2:length(x)) {
    foo <- phi %*% Gamma %*% diag(allprobs[t, ]) # phi_t-1 * Gamma * allprobs
    l <- l + log(sum(foo)) # updated log-likelihood
    phi <- foo / sum(foo) # phi_t
  }
  return(-l) # return -llk in order to maximize llk with the minimizer nlm()
}



llk(theta.star, x = data$value_kWh) # works [1] 96122.46



## Run optimizer with random starting values

# takes a while (2-3 mins)
mod <- nlm(llk, theta.star, x = data$value_kWh, print.level = 2) 
# Relative gradient close to zero. Current iterate is probably solution.

mod #  mod$minimum 4511.744
mod$estimate 
# [1]   3.8296841   1.8906948   0.9831892   2.5044483 -23.3020927  -2.1033617  -1.1639394  -0.4569191   0.7425333
# [10]  -2.6701853  -1.5727553  -0.0368120

theta.star <- mod$estimate
Gamma <- transform_3x3(theta.star[1:6]) # retransform probabilities - to be in [0, 1], each rowsum 1
# retransform mean and standard deviations - to be positive
mu <- exp(theta.star[7:9])
sigma <- exp(theta.star[10:12])





# **************** Finding optimum by running several randomly initialized models - each takes ~3 mins
llks <- rep(NA, 5) # we will run 5 models
mods <- list()
set.seed(123)
for (k in 1:5) { 
  # next line: 6 random nums ~ 0 as initial transition param, 6 random nums between 0.05 and 5 as the parameters
  theta.star <- c(rnorm(6, mean = 0, sd = 1), log(runif(6, 0.05, 5)))
  mods[[k]] <- nlm(llk, theta.star, x = data$value_kWh) 
  llks[k] <- -mods[[k]]$minimum
  print(k)
}

llks # all but one are the same. -> probably a global maximum
# [1]  -4511.753  -4511.744  -4511.744 -10552.696  -4511.744


## Backtransform parameters
theta.star <- mods[[which(llks == max(llks))]]$estimate # choose best result (highest likelihood)
Gamma <- transform_3x3(theta.star[1:6]) # retransform probabilities - to be in [0, 1], each rowsum 1
# retransform mean and standard deviations - to be positive
mu <- exp(theta.star[7:9])
sigma <- exp(theta.star[10:12])



# **************************** VISUALIZING THE RESULT
# get stationary distribution of the resulting MC
# combining the conditions of delta*Gamma=delta and sum of deltas=1 to one linear system
delta <- solve(t(diag(3) - Gamma + 1), rep(1, 3)) # (diagonal matrix - gamma matrix + 1 to every value)transposed = (1, 1, 1)transposed


## Plot results
hist(data$value_kWh, breaks = 60, col = "grey", 
     xlab = "Power usage in kWh", main = "Histogram", prob = T, 
     border=0)

x <- seq(0, 6, by = 0.01)
lines(x, dgamma(x, shape = mu[1]^2 / sigma[1]^2,
                scale = sigma[1]^2 / mu[1]) * delta[1], col = "orange", lwd = 3)
lines(x, dgamma(x, shape = mu[2]^2 / sigma[2]^2,
                scale = sigma[2]^2 / mu[2]) * delta[2], col = "red", lwd = 3)
lines(x, dgamma(x, shape = mu[3]^2 / sigma[3]^2,
                scale = sigma[3]^2 / mu[3]) * delta[3], col = "purple", lwd = 3)

lines(x, dgamma(x, shape = mu[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu[1]) * delta[1] +
        dgamma(x, shape = mu[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu[2]) * delta[2] +
        dgamma(x, shape = mu[3]^2 / sigma[3]^2, scale = sigma[3]^2 / mu[3]) * delta[3],
      col = "black", lty = 2, lwd = 3
)
legend("topright", c("state 1", "state 2", "state 3", "mixture"),
       col = c("orange", "red", "purple", "black"),
       lty = c(1, 1, 1, 2), lwd = c(3, 3, 3, 3)
)





# *******************************************************************************************
# *******************************************************************************************






### Random search with LaMa HMMs 2-5 states, 10 models for each, to find the optimum

library(LaMa)

nll <- function(theta.star){ # negative log likelihood function with LaMa
  getAll(theta.star, dat)
  
  Gamma <- tpm(logit_Gamma)
  delta <- stationary(Gamma)
  sigma <- exp(log_sigma)
  mu <- exp(log_mu)
  
  REPORT(mu); ADREPORT(mu)
  REPORT(sigma); ADREPORT(sigma)
  REPORT(Gamma)
  REPORT(delta)
  
  
  allprobs <- matrix(1, length(dat$d), dat$N)
  ind <- which(!is.na(dat$d))
  for (j in 1:dat$N){
    allprobs[ind, j] <- dgamma2(dat$d[ind],mu[j],sigma[j])  # dnorm(x[ind], mu[j], sigma[j])
  }
  
  REPORT(allprobs)
  
  -forward(delta, Gamma, allprobs)
}

# *************************************************
# function to generate starting parameters
random_par <- function(N, x){ 
  m <- mean(x)
  s <- sd(x)
  list(
    logit_Gamma = rnorm(N * (N - 1), 0, 2), # small values in range -4, +4 
    # so Gamma has reasonable range of values
    log_mu = log(runif(N, 0.1*m, 2*m)),
    log_sigma = log(runif(N, 0.1*s, 2*s))
  )
}


# random search results stored here
results <- list()

for(N in 2:5){ # for 2-5 number of states
  dat <- list(N = N, d = data$value_kWh)
  llks <- numeric(10) # a vector of length 10
  mods <- vector("list",10)
  
  for(i in 1:10){ # 10 models, 10 different lists of starting parameters
    par <- random_par(N, data$value_kWh)
    obj <- MakeADFun(nll, par)
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    mod <- report(obj)
    mods[[i]] <- list(opt = opt, mod = mod) # storing the nll result and model   # used opt instead of nlminb(...), bc we already have it stored
    llks[i] <- -mods[[i]]$opt$objective # transform negative llk to positive
    cat("States:",N," Run:",i, "likelihood:", llks[i], "\n") # print to console
  }
  
  results[[paste0("N",N)]] <- list(llks = llks,
                                   best = mods[[which.max(llks)]] # model with maximum llk
  )
}

# *************************************************
# VISUALIZATION

for (N in 2:5){
  cat(N, "states best (min nll) result:", results[[paste0("N", N)]]$best$opt$objective, "\n")
  best_5S_model = results[[paste0("N", N)]]$best$mod
  mod = best_5S_model
  
  mu <- mod$mu
  sigma <- mod$sigma
  Gamma <- mod$Gamma
  delta <- mod$delta
  
  x <- seq(0, max(data$value_kWh), length.out = 500)   # plotting range
  
  # histogram
  hist(
    data$value_kWh, breaks = 50, col = "grey",
    xlab = "Power usage in kWh",
    main = paste(N, "state model"),
    prob = TRUE, border = 0
  )
  
  cols <- rainbow(N) # colors for states
  
  # plot individual states
  for(j in 1:N){
    lines(x,
          dgamma(x, shape = mu[j]^2 / sigma[j]^2, scale = sigma[j]^2 / mu[j]) * delta[j],
          col = cols[j], lwd = 3)
  }
  
  # mixture density
  mixture <- rep(0,length(x))
  for(j in 1:N){
    mixture <- mixture +
      dgamma(x, shape = mu[j]^2 / sigma[j]^2, scale = sigma[j]^2 / mu[j]) * delta[j]
  }
  lines(x, mixture, col="black", lwd=3, lty=2)
  
  legend("topright",
         c(paste("state",1:N), "mixture"), col=c(cols,"black"),
         lty=c(rep(1,N),2), lwd=3)
}


# *************************************************
## AIC and BIC for the models:

T <- sum(!is.na(data$value_kWh))
# this will be a nice table
model_selection <- data.frame(States = 2:5, LogLik = NA,
                              NumOfParameters = NA, AIC = NA, BIC = NA)

for (N in 2:5) {
  best <- results[[paste0("N", N)]]$best
  llk <- -best$opt$objective
  # number of par k is Gamma pars, mus, sigmas, N is the number of states. k=N(N−1)+N+N=N(N+1)
  k <- N * (N + 1)
  
  model_selection[model_selection$States == N, ] <- list(
    N, llk,
    k, -2 * llk + 2 * k, -2 * llk + log(T) * k)
}

model_selection 


# *******************************************************************************************
### State-decoding 

# Global decoding - Viterbi
states <- viterbi(mod = mod)

col = c("red", "purple", "blue", "orange", "green")
mod$mu #  1.2133267 0.3430480 2.6544948 0.5965250 0.2497793
# so: green, purple, orange, red, blue would be expected

par(mfrow=c(1, 1))

idx <- 1:35900 # overall
idx <- 5137:(5136 + 8760) # 2017

idx <- 5137:(5136 + 48) # 2 days in january
idx <- (5137 + 4300):(5137 + 4300 + 48) # 2 days in the summer
plot(data$start_datetime[idx], data$value_kWh[idx], type = "n")  # empty plot with axes
for (i in seq_along(idx)) {
  points(data$start_datetime[idx[i]], data$value_kWh[idx[i]], col = col[states[idx[i]]], pch = 16)
}



# *************************************************
## RESIDUALS FOR EACH BEST MODEL
for (N in 2:5){
  best_model = results[[paste0("N", N)]]$best$mod
  mod = best_model
  
  pres <- pseudo_res(obs = data$value_kWh,
                     dist = "gamma2", # have to tell which distribution
                     par = list(mean = mod$mu, sd = mod$sigma),
                     mod = mod)
  plot(pres)
}

# autocorrelation terrible in all but gets better with more states, 
# pseudo-res normally distributed, QQ plot gets better with mode states

# that means we extracted not all info, so we can try to add covariates



# *******************************************************************************************
# *******************************************************************************************

### Covariates
# moy

nll_moy <- function(theta.star){  # negative log likelihood function + covariates with LaMa
  getAll(theta.star, dat)
  
  Gamma <- tpm_g(Z, beta); 
  Delta = stationary_p(Gamma); 
  
  delta <- Delta[1, ]
  
  sigma <- exp(log_sigma)
  mu <- exp(log_mu)
  
  REPORT(mu); ADREPORT(mu)
  REPORT(sigma); ADREPORT(sigma)
  REPORT(Gamma); ADREPORT(Gamma)
  REPORT(Delta); ADREPORT(Delta)
  
  allprobs <- matrix(1, length(dat$d), dat$N)
  ind <- which(!is.na(dat$d))
  for (j in 1:dat$N){
    allprobs[ind, j] <- dgamma2(dat$d[ind],mu[j],sigma[j])
  }
  
  REPORT(allprobs)
  
  -forward_g(delta, Gamma, allprobs)
}


# function to generate starting parameters for HMM with covariates
random_par_cov <- function(N, x){ 
  m <- mean(x)
  s <- sd(x)
  
  list(
    beta = cbind(rep(-1, N*(N - 1)), matrix(0, N*(N - 1), 2)),
    log_mu = log(runif(N, 0.1*m, 2*m)),
    log_sigma = log(runif(N, 0.1*s, 2*s))
  )
}

results_c <- list()

for(N in 2:5){ # for 2-5 number of states
  Z = cbind(
    sin(2 * pi * data$month_of_year / 12),
    cos(2 * pi * data$month_of_year / 12))  # harmonize our seasonality
  
  dat <- list(N = N, d = data$value_kWh, Z = Z)
  llks <- numeric(10)  # a vector of length 10
  mods <- vector("list", 10)
  
  for(i in 1:10){  # 10 models, 10 different lists of starting parameters
    par <- random_par_cov(N, data$value_kWh)
    obj <- MakeADFun(nll_moy, par, silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    mod <- report(obj)
    
    mods[[i]] <- list(opt = opt, mod = mod, par_list = obj$env$parList(opt$par))  # storing the nll result and model
    llks[i] <- -opt$objective  # transform negative llk to positive
    
    cat("States:",N," Run:",i, "likelihood:", llks[i], "\n") # print to console
  }
  
  results_c[[paste0("N",N)]] <- list(llks = llks,
                                   best = mods[[which.max(llks)]] # model with maximum llk
  )
}

names(results_c)


# *************************************************
## AIC and BIC for the models:

T <- sum(!is.na(data$value_kWh))
# this will be a nice table
model_selection_cov <- data.frame(States = 2:5, LogLik = NA,
                                  NumOfParameters = NA, AIC = NA, BIC = NA)

for (N in 2:5) {
  best <- results_c[[paste0("N", N)]]$best
  llk <- -best$opt$objective
  
  # Number of parameters (k) calculation:
  # ** beta: N*(N - 1) transitions, each has 3 parameters (intercept, sin, cos) -> 3*N*(N - 1)
  # ** Gamma distrib. params.: 2 parameters (mu, sigma) per state -> 2*N
  # Total k = 3*N*(N - 1) + 2*N = 3*N^2 - N = N*(3*N - 1)
  k <- N*(3*N - 1)
  
  model_selection_cov[model_selection_cov$States == N, ] <- list(
    N, llk, k,
    -2 * llk + 2 * k,
    -2 * llk + log(T) * k)
}

model_selection_cov


# *************************************************
# VISUALIZATION for N = 5 (best model)

best_model <- results_c$N5$best  # taking best model with 5 parameters
best_par_list <- best_model$par_list


Z_mini <- cbind(
  sin(2 * pi * (1:12) / 12),
  cos(2 * pi * (1:12) / 12))

dat <- list(N = 5, d = rep(1, 12), Z = Z_mini)  # mini data for visualization of the already fitted Delta

obj_N5 <- MakeADFun(nll_moy, best_par_list, silent = TRUE)

sdr_N5 <- sdreport(obj_N5)

Delta <- as.list(sdr_N5, "Est", report = TRUE)$Delta

colors <- rainbow(5)


plot(1:12, type = "n", ylim = c(0, 1), xlim = c(1, 12),
     xlab = "Month", ylab = "State probability", 
     bty = "n", xaxt = "n")
axis(1, at = 1:12)

# drawing the graphs
for(j in 1:5) {
  est <- Delta[, j]
 
  lines(1:12, est, type = "b", col = colors[j], lwd = 2, pch = 16)
}


legend("topright", legend = paste("state", 1:5), col = colors, 
       lty = 1, lwd = 2.5, bty = "n")


# *************************************************
## RESIDUALS FOR EACH BEST MODEL
for (N in 2:5){
  best_model = results_c[[paste0("N", N)]]$best$mod
  mod = best_model
  
  pres <- pseudo_res(obs = data$value_kWh,
                     dist = "gamma2", # have to tell which distribution
                     par = list(mean = mod$mu, sd = mod$sigma),
                     mod = mod)
  plot(pres, main = N)
}



# *******************************************************************************************
# *******************************************************************************************
# tod

# we can use previous nnl_moy() as nll function

results_c2 <- list()

for(N in 2:5){ # for 2-5 number of states
  Z = cbind(
    sin(2 * pi * data$hour_of_day / 24),
    cos(2 * pi * data$hour_of_day / 24)
  )  
  
  dat <- list(N = N, d = data$value_kWh, Z = Z)
  llks <- numeric(10) 
  mods <- vector("list", 10)
  
  for(i in 1:10){ 
    par <- random_par_cov(N, data$value_kWh)
    obj <- MakeADFun(nll_moy, par, silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    mod <- report(obj)
    
    mods[[i]] <- list(opt = opt, mod = mod, par_list = obj$env$parList(opt$par)) 
    llks[i] <- -opt$objective 
    
    cat("States:", N, " Run:", i, "likelihood:", llks[i], "\n") 
  }
  
  results_c2[[paste0("N", N)]] <- list(
    llks = llks,
    best = mods[[which.max(llks)]] 
  )
}

names(results_c2)


# *************************************************
## AIC and BIC for the models:

T <- sum(!is.na(data$value_kWh))
# this will be a nice table
model_selection_cov2 <- data.frame(States = 2:5, LogLik = NA,
                                  NumOfParameters = NA, AIC = NA, BIC = NA)

for (N in 2:5) {
  best <- results_c2[[paste0("N", N)]]$best
  llk <- -best$opt$objective
  k <- N * (3 * N - 1)
  
  model_selection_cov2[model_selection_cov2$States == N, ] <- list(
    N, llk, k,
    -2 * llk + 2 * k,
    -2 * llk + log(T) * k)
}

model_selection_cov2


# *************************************************
# VISUALIZATION for N = 5 (best model)

best_model <- results_c2$N5$best  # taking best model with 5 parameters
best_par_list <- best_model$par_list


Z_mini <- cbind(
  sin(2 * pi * (0:23) / 24),
  cos(2 * pi * (0:23) / 24))

dat <- list(N = 5, d = rep(1, 24), Z = Z_mini)  # for visualization of Delta which was already found

obj_N5 <- MakeADFun(nll_moy, best_par_list, silent = TRUE)

sdr_N5 <- sdreport(obj_N5)

Delta <- as.list(sdr_N5, "Est", report = TRUE)$Delta

colors <- rainbow(5)


plot(0:23, type = "n", ylim = c(0, 1), xlim = c(0, 23),
     xlab = "Hour of day", ylab = "State probability", 
     bty = "n", xaxt = "n")
axis(1, at = 0:23)

# drawing the graphs
for(j in 1:5) {
  est <- Delta[, j]
  
  lines(0:23, est, type = "b", col = colors[j], lwd = 2, pch = 16)
}


legend("topright", legend = paste("state", 1:5), col = colors, 
       lty = 1, lwd = 2.5, bty = "n")


# *************************************************
## RESIDUALS FOR EACH BEST MODEL
for (N in 2:5){
  best_model = results_c2[[paste0("N", N)]]$best$mod
  mod = best_model
  
  pres <- pseudo_res(obs = data$value_kWh,
                     dist = "gamma2", # have to tell which distribution
                     par = list(mean = mod$mu, sd = mod$sigma),
                     mod = mod)
  plot(pres, main = N)
}
