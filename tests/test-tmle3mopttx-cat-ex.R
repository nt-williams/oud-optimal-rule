library(lmtp)
library(sl3)

box::use(../R/blip[estimate_blip_multi_sl3])

progressr::handlers(global = TRUE)

data("data_cat_realistic", package = "tmle3mopttx")

a <- "A"
y <- "Y"
w <- paste0("W", 1:4)

data_cat_realistic[, A := factor(A)]

Q_learner <- make_learner_stack(
  Lrnr_mean,
  Lrnr_glm_fast,
  list(Lrnr_xgboost, nrounds = 50),
  list(Lrnr_xgboost, nrounds = 100),
  list(Lrnr_xgboost, nrounds = 300)
)

g_learner <- make_learner_stack(
  Lrnr_mean,
  list(Lrnr_xgboost, nrounds = 100),
  list(Lrnr_xgboost, nrounds = 300)
)

mv_learners <- lapply(
  list(Lrnr_xgboost$new(nrounds = 50),
       Lrnr_xgboost$new(nrounds = 100),
       Lrnr_xgboost$new(nrounds = 300),
       Lrnr_mean$new(),
       Lrnr_glm_fast$new()),
  \(x) sl3::make_learner(sl3::Lrnr_multivariate, x)
)

b_learner <- sl3::make_learner(sl3::Stack, mv_learners)

policy_factory <- function(x) {
  \(data, trt) factor(rep(x, nrow(data)), levels = 1:3)
}

onestep <- purrr::partial(
  lmtp_sdr,
  data = data_cat_realistic,
  trt = a,
  outcome = y,
  baseline = w,
  folds = 10,
  .SL_folds = 10,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)

# obtain the EIFs for the counterfactual TSM under
# different trts using one-step estimator
tsm_1 <- onestep(shift = policy_factory(1))
tsm_2 <- onestep(shift = policy_factory(2))
tsm_3 <- onestep(shift = policy_factory(3))

type_2_blip <- function(...) {
  objs <- list(...)
  purrr::map_dfc(objs, \(x) x$eif) |>
    setNames(paste0("blip", seq_along(objs))) |>
    as.data.frame() |>
    (\(x) x - rowMeans(x))()
}

# estimate the blip function
data_cat_realistic <- cbind(data_cat_realistic, type_2_blip(tsm_1, tsm_2, tsm_3))
blip_estimates <- estimate_blip_multi_sl3(data_cat_realistic, w, paste0("blip", 1:3), b_learner, 10)
data_cat_dv <- data.table::copy(data_cat_realistic)

# assign trt based on optimal rule
rule <-
  c(1, 2, 3)[max.col(blip_estimates)] |>
  factor(levels = 1:3)

data_cat_dv[, A := rule]
data_cat_dv[, A := factor(tmle_spec$return_rule, levels = 1:3)]

# estimate the TSM under the optimal rule
optimal <- lmtp_tmle(
  data_cat_realistic, a, y, w,
  shifted = data_cat_dv,
  folds = 10,
  .SL_folds = 10,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)
