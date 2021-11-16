library(lmtp)
library(future)

box::use(../R/blip, ../R/utils[...], here[here], glue[glue])

model <- "super-learner-type1"

progressr::handlers(global = TRUE)

# importing imputed data and "CATEs"
oud <- readRDS(here("data", "drv", "imputed-no-27bup.rds"))

glue("estimated-{model}-rules-no27bup.rds") |>
  (\(x) here("data", "drv", x))() |>
  readRDS() -> rules

oud_optimal_rule <- map2(oud, rules, \(x, y) dplyr::mutate(x, {{ a }} := y))

plan(multisession)

# estimating tsm for 12-week relapse under the optimal rule
optimal_tsms <- map2(
  oud, oud_optimal_rule,
  function(.x, .y) {
    lmtp_sdr(
      .x, a, y, w,
      shifted = .y,
      folds = 10,
      .SL_folds = 10,
      learners_outcome = Q_learner,
      learners_trt = g_learner
    )
  }
)

plan(sequential)

saveRDS(optimal_tsms, here("data", "drv", glue("optimal-tsms-{model}-no27bup.rds")))
