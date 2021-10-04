lmtp_type <- function(x) {
    structure(x, class = class(x[[1]]))
}

rubins_rules <- \(x, ...) UseMethod("rubins_rules", lmtp_type(x))

rubins_rules.lmtp <- function(lmtps, alpha = 0.05) {
    thetas <- unlist(lapply(lmtps, \(x) x$theta))

    vw <- mean(unlist(lapply(lmtps, \(x) x$standard_error^2)))
    vb <- stats::var(thetas)

    theta <- mean(thetas)
    se <- pooled_se(vw, vb, length(thetas))

    list(
        theta = theta,
        se = se,
        alpha = alpha,
        ci = theta + (stats::qnorm(alpha / 2) * c(1, -1)) * se
    )
}

rubins_rules.lmtp_contrast <- function(lmtps, alpha = 0.05) {
    thetas <- unlist(lapply(lmtps, \(x) x$vals$theta))

    vw <- mean(unlist(lapply(lmtps, \(x) x$vals$std.error^2)))
    vb <- stats::var(thetas)

    theta <- mean(thetas)
    se <- pooled_se(vw, vb, length(thetas))

    list(
        theta = theta,
        se = se,
        alpha = alpha,
        ci = theta + (stats::qnorm(alpha / 2) * c(1, -1)) * se
    )
}

pooled_se <- \(vw, vb, m) sqrt(vw + vb + (vb / m))
