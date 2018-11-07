# Setup and read data ---------------
library(hydroGOF)   # For calculating performance scores
library(magrittr)   # To use the pipe operator
library(data.table) # For data wrangling
library(ggplot2)    # For plotting
library(cowplot)    # For plotting
library(doParallel) # For parallel computing
source('postprocess_functions.R')
source('LDS_reconstruction.R')

Qa <- fread('Data/Ping_annual.csv')
mu <- mean(log(Qa$Qa))
y <- t(c(rep(NA, 321), log(Qa$Qa)-mu))
U <- t(read.csv('Data/Ping_PCs.csv', header = F))[c(1,3,4,6,9,11,12), ] # Selected PCs

# Gridded starting conditions
grid <- expand.grid(seq(0.1, 0.9, 0.1), 
                    seq(-1, 1, 0.1),
                    seq(0.1, 1, 0.1),
                    seq(-1, 1, 0.1)) %>% 
    as.matrix() %>% 
    split(.,1:nrow(.), drop = TRUE)

# Setup parallel computing
nbCores <- detectCores()
cl <- makeCluster(nbCores)
registerDoParallel(cl)
invisible(clusterCall(cl, function() { source('EM_perfect_prediction.R')}))
invisible(clusterCall(cl, function() { source('LDS_reconstruction.R') }))

# LDS_V0: PERFECT PREDICTION ----------------------------

source('EM_perfect_prediction.R')

# LDS with 'proper' state transition ---------------------

u <- v <- U

# Learn the model with full data
# gridded starting conditions
system.time(
    lds <- learnLDS_restart(y, u, v, mu, grid, return.init = TRUE, parallel = TRUE)   
    # max lik: 1136.482
    # init: 0.9, 1, 1, 1
)
# randomized starting conditions
set.seed(100)
system.time(
    lds <- learnLDS_restart(y, u, v, mu, num.restart = 1000, return.init = TRUE, parallel = TRUE)  
    # max lik: 1039.483
    # init: 0.8353482 0.9070261 0.9295341 0.8448843
)
plot(lds$model$fit$X[1,], type = 'l')

# Test several different initial conditions
fit <- learnLDS(y, u, v, mu, c(0.4, -0.03, 0.2, -0.03))
fit$theta
fit$lik #438.4051
plot(fit$fit$X[1, ], type = 'l')

fit <- learnLDS(y, u, v, mu, c(0.9, 0.5, 0.4, 1))
fit$theta
fit$lik #438.4492
plot(fit$fit$X[1, ], type = 'l')

fit <- learnLDS(y, u, v, mu, c(0.9, -1, 1, -1))
fit$theta
fit$lik #1006.967
plot(fit$fit$X[1, ], type = 'l')

fit <- learnLDS(y, u, v, mu, c(0.9, 1, 1, 1))
fit$theta
fit$lik #1136.482
plot(fit$fit$X[1, ], type = 'l')

fit <- learnLDS(y, u, v, mu, c(0.9, 1, 1, 1))
fit$theta
fit$lik #1136.482
plot(fit$fit$X[1, ], type = 'l')

fit <- learnLDS(y, u, v, mu, c(0.84, 0.91, 0.93, 0.84)) # We got this from 1000 randomized run with seed 100
fit$theta
fit$lik #1042.939
plot(fit$fit$X[1, ], type = 'l')

# LDS with 'improper' state transtiion --------------------------------

u <- cbind(U[, 2:406], rep(0, nrow(U)))
v <- U

# Learn the model with full data
# gridded starting conditions
system.time(
    lds <- learnLDS_restart(y, u, v, mu, grid, return.init = TRUE, parallel = TRUE)   
    # max lik: 
    # init: 
)

# randomized starting conditions
set.seed(100)
system.time(
    lds <- learnLDS_restart(y, u, v, mu, num.restart = 1000, return.init = TRUE, parallel = TRUE)  
    # max lik: 900.2447
    # init: 0.1074170 -0.9862498  0.6861162  0.7610216
)
plot(lds$model$fit$X[1,], type = 'l')

# Learn this model
lds <- learnLDS(y, u, v, mu, init = c(0.11, -0.99, 0.69, 0.76))  
plot(lds$fit$X[1,], type = 'l')

# LDS_V1: CORRECT DERIVATION -------------------------

source('EM.R')
u <- v <- U

# randomized starting conditions
set.seed(100)
system.time(
    lds <- learnLDS_restart(y, u, v, mu, num.restart = 1000, return.init = TRUE, parallel = TRUE)  
    # max lik: 
    # init: 
)
plot(lds$model$fit$X[1,], type = 'l')

# --------------------------
# Plots
plot_inst_period(lds)
plot_full_period(lds)
# Benchmark
source('benchmark.R')
pc <- fread('Data/Ping_PCs.csv')
bm <- PCR_reconstruction(pc, target, 1600:2005)
bm$metrics

dt1 <- lds$metrics.dist
dt1$SNR <- NULL
dt1$model <- 'LDS'

dt2 <- bm$metrics.dist
dt2$model <- 'PCR'

dt <- rbind(dt1, dt2)
dt <- melt(dt, id.vars = 'model', variable.name = 'metric')

ggplot(dt, aes(model, value)) +
    geom_jitter(aes(group = model), colour = 'gray') +
    stat_summary(geom = 'point', colour = 'red', fun.y = mean) +
    facet_wrap(~metric, scale = 'free')
