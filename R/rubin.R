rubins_rules <- function(lmtps, alpha = 0.05) {
    thetas <- unlist(lapply(lmtps, \(x) x$theta))

    vw <- mean(unlist(lapply(lmtps, \(x) x$standard_error^2)))
    vb <- stats::var(thetas)

    theta <- mean(thetas)
    se <- pooled_se(vw, vb, length(thetas))

    list(
        theta = theta,
        se = se,
        alpha = alpha,
        ci = theta + c(stats::qnorm(0.05 / 2),
                       stats::qnorm(0.05 / 2, lower.tail = FALSE))*se
    )
}

pooled_se <- \(vw, vb, m) sqrt(vw + vb + (vb / m))
