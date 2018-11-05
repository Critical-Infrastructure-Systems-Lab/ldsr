#' Calculate some metrics from reconstructed and observed streamflow
#'
#' @param Qa.hat reconstructed
#' @param Qa observed
#' @param z indices of left-out values in cross validation
#' @export
calculate_metrics <- function(Qa.hat, Qa, z) {

    data.table(
        R2    = 1 - var(Qa.hat[-z] - Qa[-z], na.rm = TRUE) / var(Qa[-z], na.rm = TRUE),
        RE    = 1 - ssq(Qa.hat[z], Qa[z]) / ssq(rep(mean(Qa[-z], na.rm = TRUE), length(z)), Qa[z]),
        CE    = NSE(Qa.hat[z], Qa[z]),
        nRMSE = rmse(Qa.hat[z], Qa[z]) / mean(Qa, na.rm = T),
        KGE   = KGE(Qa.hat[z], Qa[z])
    )
}

#' Plot reconstruction results in the instrumental period or full period
#'
#' @param s Reconstruction results for a site
#' @param target Reconstruction target
#' @export
plot_reconstruction <- function(s, target, period) {

    if (period == 'inst') {
        s <- s[year %in% target$year]
    }
    p <- ggplot(s) +
        geom_ribbon(aes(year, ymin = Xl, ymax = Xu), fill = 'gray90') +
        geom_line(aes(year, X)) +
        labs(x = NULL, y = 'Flow regime [-]') +
        theme(legend.position = 'none',
              legend.key.width = unit(2, 'cm'))

    q <- ggplot(s) +
        geom_ribbon(aes(year, ymin = Ql, ymax = Qu), fill = 'gray90') +
        geom_line(aes(year, Q, colour = 'LDS')) +
        geom_line(aes(year, Qa, colour = 'Instrumental'), data = target) +
        labs(x = NULL, y = 'Annual streamflow [million m³]') +
        theme(legend.position = 'none',
              legend.key.width = unit(2, 'cm'))

    plot_grid(p, q, align = 'v', ncol = 1, rel_heights = c(1,1))
}

plot_cv <- function(cv.result, case.name) {

    metrics.dist <- cv.result$metrics.dist
    metrics.dist2 <- cv.result$metrics.dist2
    all.Q <- cv.result$all.Q

    # TODO: fix this, maybe copy from another folder
    #cvQ <-

    metrics.dist[, rep := 1:.N]
    metrics.dist2[, rep := 1:.N]

    p <- ggplot(all.Q) +
        geom_line(aes(year, Q, colour = 'LDS'), size = 0.1) +
        geom_line(aes(year, Qa, colour = 'Inst'), data = Qa, size = 0.1) +
        geom_point(aes(year, cvQ, colour = 'LDS'), data = cvQ, size = 0.1) +
        geom_point(aes(year, cvObs, colour = 'Inst'), data = cvQ, size = 0.1) +
        geom_text(aes(x = 1922, y = 4000, label = paste0('CE = ', round(CE, 2))),
                  data = metrics.dist,
                  hjust = 'left') +
        facet_wrap(~rep)
    ggsave(paste0('cv_results_', case.name, 'with_update.pdf'),
           p,
           width = 420, height = 297, units = 'mm')

    q <- ggplot(all.Q) +
        geom_line(aes(year, Q2, colour = 'LDS'), size = 0.1) +
        geom_line(aes(year, Qa, colour = 'Inst'), data = Qa, size = 0.1) +
        geom_point(aes(year, cvQ2, colour = 'LDS'), data = cvQ, size = 0.1) +
        geom_point(aes(year, cvObs, colour = 'Inst'), data = cvQ, size = 0.1) +
        geom_text(aes(x = 1922, y = 4000, label = paste0('CE = ', round(CE, 2))),
                  data = metrics.dist,
                  hjust = 'left') +
        facet_wrap(~rep)
    ggsave(paste0('cv_results_', case.name, 'without_update.pdf'),
           q,
           width = 420, height = 297, units = 'mm')
}
