library(glmnet)

blips <- readRDS(here::here("data", "drv", "estimated-adaptLASSO-type1-blips.rds"))

selected <- function(x, blip, prop = TRUE, name) {
  out <-
    lapply(x, \(x) lapply(x$selected, \(y) y[[blip]])) |>
    unlist() |>
    lapply(\(x) as.matrix(x)) |>
    (\(x) if (prop) lapply(x, \(y) y != 0) else x)() |>
    (\(x) Reduce(`+`, x) / 50)()

  colnames(out) <- name
  out
}

prop_selected <- cbind(
  selected(blips, "blip1", name = "Methadone vs. Bupenorphine"),
  selected(blips, "blip2", name = "Naltrexone vs. Bupenorpine")
)[-1, , drop = FALSE]

mean_coef <- cbind(
  selected(blips, "blip1", FALSE, name = "Methadone vs. Bupenorphine"),
  selected(blips, "blip2", FALSE, name = "Naltrexone vs. Bupenorpine")
) |>
  round(4)
