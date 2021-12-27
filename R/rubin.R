lmtp_type <- \(x) structure(x, class = class(x[[1]]))

rubins_rules <- \(x, ...) UseMethod("rubins_rules", lmtp_type(x))

rubins_rules.numeric <- function(eifs, label, alpha = 0.05) {
    thetas <- unlist(lapply(eifs, mean))
    vw <- mean(unlist(lapply(eifs, \(x) stats::var(x) / length(x))))
    vb <- stats::var(thetas)
    theta <- mean(thetas)
    se <- pooled_se(vw, vb, length(thetas))
    data.frame(
        label = label,
        theta = theta,
        se = se,
        alpha = alpha,
        conf.low = theta + stats::qnorm(alpha / 2) * se,
        conf.high = theta + (stats::qnorm(alpha / 2) * -1) * se
    )
}

rubins_rules.lmtp <- function(lmtps, label, alpha = 0.05) {
    thetas <- unlist(lapply(lmtps, \(x) x$theta))

    vw <- mean(unlist(lapply(lmtps, \(x) x$standard_error^2)))
    vb <- stats::var(thetas)

    theta <- mean(thetas)
    se <- pooled_se(vw, vb, length(thetas))

    data.frame(
        label = label,
        theta = theta,
        se = se,
        alpha = alpha,
        conf.low = theta + stats::qnorm(alpha / 2) * se,
        conf.high = theta + (stats::qnorm(alpha / 2) * -1) * se
    )
}

rubins_rules.lmtp_contrast <- function(lmtps, label, alpha = 0.05) {
    thetas <- unlist(lapply(lmtps, \(x) x$vals$theta))

    vw <- mean(unlist(lapply(lmtps, \(x) x$vals$std.error^2)))
    vb <- stats::var(thetas)

    theta <- mean(thetas)
    se <- pooled_se(vw, vb, length(thetas))

    data.frame(
        label = label,
        theta = theta,
        se = se,
        alpha = alpha,
        conf.low = theta + stats::qnorm(alpha / 2) * se,
        conf.high = theta + (stats::qnorm(alpha / 2) * -1) * se
    )
}

pooled_se <- \(vw, vb, m) sqrt(vw + vb + (vb / m))
