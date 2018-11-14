bm <- PCR_reconstruction(P1annual, P1pc)
bm$metrics
bm2 <- PCR_reconstruction(P1annual, P1pc[,c(1,4,6,9,11,12,13)])
bm2$metrics
