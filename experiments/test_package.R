library(ldsr)
# Check that package is self-contained so I shouldn't need to load data.table,
# doParallel, ggplot2 and magrittr to run package functions

Qa <- Ping.annual
u <- Ping.PC
v <- Ping.PC

# Sequential learning
# With init
fit <- LDS_reconstruction(Qa, u, v,
                          init = c(0.4, rep(-0.03, 7), 0.2, rep(-0.03, 7)),
                          parallel = FALSE)
# Randomized
fit <- LDS_reconstruction(Qa, u, v, num.restart = 5, parallel = FALSE)
# Parallel learning
# With init
fit <- LDS_reconstruction(Qa, u, v,
                          init = c(0.4, rep(-0.03, 7), 0.2, rep(-0.03, 7)),
                          parallel = TRUE)
# Randomized
fit <- LDS_reconstruction(Qa, u, v, num.restart = 12, parallel = TRUE)

# Plot results
plot_reconstruction(fit$rec, Qa, 'inst')
plot_reconstruction(fit$rec, Qa, 'full')

# Sequential CV
# With init
cv <- cvLDS(Qa, u, v,
            init = c(0.4, rep(-0.03, 7), 0.2, rep(-0.03, 7)),
            n.reps = 5)
# Ranomized restarts
cv <- cvLDS(Qa, u, v, num.restart = 5, n.reps = 5)
# Parallel CV
# With init
cv <- cvLDS(Qa, u, v,
            init = c(0.4, rep(-0.03, 7), 0.2, rep(-0.03, 7)),
            n.reps = 12, parallel = TRUE)
# Randomized
cv <- cvLDS(Qa, u, v, num.restart = 8, n.reps = 5, parallel = TRUE)
