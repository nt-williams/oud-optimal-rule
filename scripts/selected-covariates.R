library(glmnet)

blips <- readRDS(here::here("data", "drv", "estimated-lm-type1-blips.rds"))

selected <- function(x, blip, prop = TRUE, name) {
  out <-
    lapply(x, \(x) lapply(x$coef, \(y) y[, blip, drop = FALSE])) |>
    unlist(recursive = FALSE) |>
    lapply(\(x) as.matrix(x)) |>
    (\(x) if (prop) lapply(x, \(y) y != 0) else x)() |>
    (\(x) Reduce(`+`, x) / 50)()

  colnames(out) <- name
  out
}

prop_selected <- cbind(
  selected(blips, "blip1", name = "Methadone vs. Buprenorphine"),
  selected(blips, "blip2", name = "Naltrexone vs. Buprenorpine")
)[-1, , drop = FALSE]

mean_coef <- cbind(
  selected(blips, "blip1", FALSE, name = "Methadone vs. Buprenorphine"),
  selected(blips, "blip2", FALSE, name = "Naltrexone vs. Buprenorpine")
) |>
  round(4)
