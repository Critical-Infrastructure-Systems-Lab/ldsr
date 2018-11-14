library(ldsr)
library(data.table)
library(dplyr)
u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])
set.seed(100)
lds1 <- LDS_reconstruction(P1annual, u, v, num.restarts = 1000)
saveRDS(lds1, 'lds1000.RDS')
set.seed(100)
lds2 <- LDS_reconstruction(P1annual, u, v, num.restarts = 10000)
saveRDS(lds2, 'lds10000.RDS')

lds1$lik
lds2$lik

dt <- rbind(lds1$rec %>% mutate(case = '1'),
            lds2$rec %>% mutate(case = '2'))

ggplot(dt) +
    geom_line(aes(year, X, colour = case))

ggplot(dt) +
    geom_line(aes(year, Q, colour = case))

plot_reconstruction(lds1$rec, P1annual, 'inst')
plot_reconstruction(lds2$rec, P1annual, 'inst')

# ---------------------------------------
set.seed(200)
Z <- replicate(100, sample(1:85, 9), simplify = F)
# Build and cross-validate a PCR model
bm <- PCR_reconstruction(P1annual, P1pc)
# Model performance
bm$metrics
set.seed(100)
system.time(
    cv <- cvLDS(P1annual, u, v, Z = Z, num.restarts = 1000)
)

dt1 <- bm$metrics.dist
dt1[, model := 'PCR']
dt2 <- cv$metrics.dist
dt2[, model := 'LDS']
dt <- rbind(dt1, dt2)
dt <- melt(dt, id.vars = 'model', measure.vars = 1:5, variable.name = 'metric')

ggplot(dt, aes(model, value)) +
    geom_jitter(colour = 'gray70', width = 0.2) +
    stat_summary(geom = 'point', fun.y = 'mean', colour = 'red') +
    facet_wrap(~metric, scales = 'free') +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = 'gray70', fill = NA)
    )
