# Setup -----------------------
library(ldsr)
# Check that package is self-contained so I shouldn't need to load data.table,
# doParallel, ggplot2 and magrittr to run package functions

Qa <- Ping.annual
u <- Ping.PC
v <- Ping.PC

fit <- LDS_reconstruction(Qa, u, v, num.restart = 2, parallel = TRUE)

p <- plot_reconstruction(fit$rec, Qa, 'full')
print(p)
cowplot::ggsave('Test.pdf', p)
cv <- cvLDS(Qa, u, v, init = c(0.4, -0.03, 0.2, -0.03), n.reps = 2)
