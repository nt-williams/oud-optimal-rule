library(glmnet)

blips <- readRDS(here::here("data", "drv", "estimated-adaptLASSO-type1-blips.rds"))

blip1 <-
  lapply(blips, \(x) lapply(x$selected, \(y) y$blip1)) |>
  unlist() |>
  lapply(\(x) as.matrix(x)) |>
  lapply(\(x) x != 0) |>
  (\(x) Reduce(`+`, x) / 50)()

lapply(blips, \(x) lapply(x$selected, \(y) y$blip1)) |>
  unlist() |>
  lapply(\(x) as.matrix(x)) |>
  (\(x) Reduce(`+`, x) / 50)() |>
  round(4)

blip2 <-
  lapply(blips, \(x) lapply(x$selected, \(y) y$blip2)) |>
  unlist() |>
  lapply(\(x) as.matrix(x)) |>
  lapply(\(x) x != 0) |>
  (\(x) Reduce(`+`, x) / 50)()

lapply(blips, \(x) lapply(x$selected, \(y) y$blip2)) |>
  unlist() |>
  lapply(\(x) as.matrix(x)) |>
  (\(x) Reduce(`+`, x) / 50)() |>
  round(4)
