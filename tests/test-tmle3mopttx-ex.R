library(lmtp)
library(sl3)

box::use(../R/blip[estimate_blip])

progressr::handlers(global = TRUE)

data("data_bin", package = "tmle3mopttx")

a <- "A"
y <- "Y"
w <- paste0("W", 1:3)

Q_learner <- make_learner_stack("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_xgboost")
g_learner <- make_learner_stack("Lrnr_mean", "Lrnr_glm_fast", "Lrnr_xgboost")
b_learner <- make_learner_stack("Lrnr_glm_fast", "Lrnr_xgboost")

# obtain the EIFs for the counterfactual TSM under different trts using one-step estimator
tsm_1 <- lmtp_sdr(
  data_bin, a, y, w,
  shift = static_binary_on,
  folds = 5,
  .SL_folds = 5,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)

tsm_0 <- lmtp_sdr(
  data_bin, a, y, w,
  shift = static_binary_off,
  folds = 5,
  .SL_folds = 5,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)

# estimate the blip function
data_bin$blip <- tsm_1$eif - tsm_0$eif
blip_estimates <- estimate_blip(data_bin, w, "blip", b_learner, 10)

data_bin_dv <- data.table::copy(data_bin)

# assign trt based on optimal rule
data_bin_dv[, A := as.numeric(blip_estimates > 0)]

# estimate the TSM under the optimal rule
optimal <- lmtp_tmle(
  data_bin, a, y, w,
  shifted = data_bin_dv,
  folds = 10,
  .SL_folds = 10,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)
