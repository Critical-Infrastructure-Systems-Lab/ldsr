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
#' @param period Either 'inst' or 'full'
#' @param benchmark A benchmark reconstruction results following the format of `PCR_reconstruction`
#' @export
plot_reconstruction <- function(s, target, period, benchmark = NULL) {

    if (period == 'inst') s <- s[year %in% target$year]
    p <- ggplot(s) +
        theme_cowplot() +
        panel_border('black') +
        geom_ribbon(aes(year, ymin = Ql, ymax = Qu), alpha = 0.2) +
        geom_line(aes(year, Q, colour = 'LDS', linetype = 'LDS'), size = 0.4) +
        labs(x = 'Year', y = 'Annual streamflow, million m\u00B3')

    if (period == 'inst') {
        p <- p + geom_line(aes(year, Qa, colour = 'Instrumental', linetype = 'Instrumental'),
                           data = target, size = 0.4)

        if (!is.null(benchmark)) {
            benchmark <- benchmark[year %in% target$year]
            p <- p +
                geom_line(aes(year, Q, colour = 'Benchmark', linetype = 'Benchmark'),
                          data = benchmark) +
                scale_colour_manual(name = "Model", values = c('black', 'darkorange', 'black')) +
                scale_linetype_manual(name = "Model", values = c(2, 1, 1))
        } else
            p <- p +
                scale_colour_manual(name = "Model", values = c('darkorange', 'black')) +
                scale_linetype_manual(name = "Model", values = c(1, 1))

        p <- p + theme(legend.position = 'top',
                       legend.title = element_blank())
    } else {# Full period
        if (!is.null(benchmark))
            p <- p + geom_line(aes(year, Q, colour = 'Benchmark', linetype = 'Benchmark'),
                               data = benchmark) +
                scale_colour_manual(name = "Model", values = c('blue', 'black')) +
                scale_linetype_manual(name = "Model", values = c(2, 1)) +
                theme(legend.position = 'top',
                      legend.title = element_blank())
        else
            p <- p + theme(legend.position = 'none')
    }

    q <- ggplot(s) +
        geom_hline(yintercept = 0, colour = 'black', size = 0.2) +
        geom_ribbon(aes(year, ymin = Xl, ymax = Xu), alpha = 0.2) +
        geom_line(aes(year, X), size = 0.4) +
        labs(x = 'Year', y = 'Flow regime state') +
        theme_cowplot() +
        panel_border('black')

    plot_grid(p, q, nrow = 2,
              rel_heights = {if (p$theme$legend.position == 'none') c(1,1) else c(1, 0.75)},
              labels = c('(a)', '(b)'))
}

# plot_cv <- function(cv.result, case.name) {
#
#     metrics.dist <- cv.result$metrics.dist
#     metrics.dist2 <- cv.result$metrics.dist2
#     all.Q <- cv.result$all.Q
#
#     # TODO: fix this, maybe copy from another folder
#     #cvQ <-
#
#     metrics.dist[, rep := 1:.N]
#     metrics.dist2[, rep := 1:.N]
#
#     p <- ggplot(all.Q) +
#         geom_line(aes(year, Q, colour = 'LDS'), size = 0.1) +
#         geom_line(aes(year, Qa, colour = 'Inst'), data = Qa, size = 0.1) +
#         geom_point(aes(year, cvQ, colour = 'LDS'), data = cvQ, size = 0.1) +
#         geom_point(aes(year, cvObs, colour = 'Inst'), data = cvQ, size = 0.1) +
#         geom_text(aes(x = 1922, y = 4000, label = paste0('CE = ', round(CE, 2))),
#                   data = metrics.dist,
#                   hjust = 'left') +
#         facet_wrap(~rep)
#     ggsave(paste0('cv_results_', case.name, 'with_update.pdf'),
#            p,
#            width = 420, height = 297, units = 'mm')
#
#     q <- ggplot(all.Q) +
#         geom_line(aes(year, Q2, colour = 'LDS'), size = 0.1) +
#         geom_line(aes(year, Qa, colour = 'Inst'), data = Qa, size = 0.1) +
#         geom_point(aes(year, cvQ2, colour = 'LDS'), data = cvQ, size = 0.1) +
#         geom_point(aes(year, cvObs, colour = 'Inst'), data = cvQ, size = 0.1) +
#         geom_text(aes(x = 1922, y = 4000, label = paste0('CE = ', round(CE, 2))),
#                   data = metrics.dist,
#                   hjust = 'left') +
#         facet_wrap(~rep)
#     ggsave(paste0('cv_results_', case.name, 'without_update.pdf'),
#            q,
#            width = 420, height = 297, units = 'mm')
# }
