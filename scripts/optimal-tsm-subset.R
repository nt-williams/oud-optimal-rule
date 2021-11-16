library(lmtp)
library(here)

model <- "adaptLASSO-type1"

here("data/drv", glue("optimal-tsms-{model}-no27bup.rds")) |>
  readRDS() -> fits

here("data/drv", "imputed-no-27bup.rds") |>
  readRDS() |>
  (\(x) x[[1]])() |>
  (\(x) x$project == "27")() -> in27

list(
  in27_eif = lapply(fits, \(fit) fit$eif[in27]),
  no27_eif = lapply(fits, \(fit) fit$eif[!in27])
) |>
  saveRDS(here("data/drv", glue("optimal-tsms-{model}-no27bup-projectSubsets.rds")))
