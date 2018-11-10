library(ldsr)
# Check that package is self-contained so I shouldn't need to load data.table,
# doParallel, ggplot2 and magrittr to run package functions

Qa <- Ping.annual
u <- Ping.PC
v <- Ping.PC

# LDS EM ------------------------------------------------
# Sequential learning
# With init
fit <- LDS_reconstruction(Qa, u, v, method = 'EM',
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

# With reupdate
cv <- cvLDS(Qa, u, v, method = 'EM', num.restart = 8, n.reps = 5, parallel = TRUE, reupdate = TRUE)

# LDS GA ---------------------------------------------------

fit <- LDS_reconstruction(Qa, u, v, method = 'GA', niter = 10000, parallel = 8)
fit$theta
fit$lik
plot_reconstruction(fit$rec, Qa, 'full')
plot_reconstruction(fit$rec, Qa, 'inst')

cv <- cvLDS(Qa, u, v, lambda = 1, method = 'GA', num.restart = 2, n.reps = 8, parallel = TRUE)

# Package calls
fit <- LDS_reconstruction(Qa, u, v, method = 'EM')
fit <- LDS_reconstruction(Qa, u, v, method = 'GA')
fit <- LDS_reconstruction(Qa, u, v, method = 'BFGS')
fit <- LDS_reconstruction(Qa, u, v, method = 'abcd')

check <- function(x) {
    if (is.list(x))
        "list"
    else
        "not list"

}
check(list(a = 1, b = 2, c = 3))
check(c(a = 1, b = 2, c = 3))

