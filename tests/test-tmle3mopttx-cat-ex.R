library(lmtp)
library(sl3)
library(future)

box::use(../R/blip[...])

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
  folds = 3,
  .SL_folds = 3,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)

plan(multisession)

# obtain the EIFs for the counterfactual TSM under
# different trts using one-step estimator
tsm_1 <- onestep(shift = policy_factory(1))
tsm_2 <- onestep(shift = policy_factory(2))
tsm_3 <- onestep(shift = policy_factory(3))

# estimate the blip function
data_cat_realistic <- cbind(data_cat_realistic, type_1_blip(tsm_1, tsm_2, ref = tsm_3))
blip_estimates <- estimate_blip_multi_sl3(data_cat_realistic, w, paste0("blip", 1:2), b_learner, 3)

# assign trt based on optimal rule
rule <- find_optimal_rule(1:2, blip_estimates) |>
  factor(levels = 1:3)

data_cat_dv <- data.table::copy(data_cat_realistic)
data_cat_dv[, A := rule]

# estimate the TSM under the optimal rule
optimal <- lmtp_sdr(
  data_cat_realistic, a, y, w,
  shifted = data_cat_dv,
  folds = 10,
  .SL_folds = 10,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)
# LMTP Estimator: SDR
#    Trt. Policy:
#
# Population intervention effect
#       Estimate: 0.6674
#     Std. error: 0.0615
#         95% CI: (0.547, 0.7879)
