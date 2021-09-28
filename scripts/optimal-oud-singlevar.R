library(lmtp)
library(sl3)

# load functions
box::use(../R/blip[...], ../R/utils[...])

# setting up parallel processing
plan(multisession)

progressr::handlers(global = TRUE)

oud <- readRDS(here::here("data", "drv", "onestep-tsm-imputed.rds"))

# SL stack for outcome regression
Q_learner <- make_learner_stack(
  Lrnr_mean,
  Lrnr_glm_fast,
  Lrnr_earth,
  list(Lrnr_xgboost, nrounds = 50),
  list(Lrnr_xgboost, nrounds = 100),
  list(Lrnr_xgboost, nrounds = 300)
)

# SL stack for propensity score
g_learner <- make_learner_stack(
  Lrnr_mean,
  Lrnr_glm_fast,
  Lrnr_earth,
  list(Lrnr_xgboost, nrounds = 50),
  list(Lrnr_xgboost, nrounds = 100),
  list(Lrnr_xgboost, nrounds = 300)
)

ans <- list()
for (i in seq_along(oud)[1]) {
  .data <- oud[[i]]$data

  # estimate the blip function
  blips <- estimate_blip_multi_singlevar(.data, w, paste0("blip", 1:3), 10)

  # assign medicine based on optimal rule
  rule <- find_optimal_rule(1:3, blips$preds, TRUE) |> factor(levels = 1:3)

  oud_optimal <- .data
  oud_optimal[, a] <- rule

  # estimate the TSM under the optimal rule
  opt <- lmtp_sdr(
    .data, a, y, w,
    shifted = oud_optimal,
    folds = 10,
    .SL_folds = 10,
    learners_outcome = Q_learner,
    learners_trt = g_learner
  )

  ans[[i]] <- list(
    blip_estimates = blips$preds,
    selected = blips$selected,
    rule = rule,
    optimal = opt
  )
}
