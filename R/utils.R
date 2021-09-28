a <- "medicine"                                           # medicine variable
y <- "week_12_relapse"                                    # outcome variable

w <- c(                                                   # potential effect modifiers
  "sex", "age", "xrace", "hwithdraw", "alcdisorder",
  "cocdisorder", "hasBrainDamage", "hasEpilepsy",
  "hasSchiz", "hasBipolar", "hasAnxPan", "hasMajorDep",
  "bamphetamine30_base", "bcannabis30_base",
  "bbenzo30_base", "ivdrug"
)

#' Medicine shift function
drug_assignment <- function(x) {
  \(data, trt) factor(rep(x, nrow(data)), levels = c("bup", "met", "nal"))
}

#' SL stack for outcome regression
Q_learner <- sl3::make_learner_stack(
  sl3::Lrnr_mean,
  sl3::Lrnr_glm_fast,
  sl3::Lrnr_earth,
  list(sl3::Lrnr_xgboost, nrounds = 50),
  list(sl3::Lrnr_xgboost, nrounds = 100),
  list(sl3::Lrnr_xgboost, nrounds = 300)
)

#' SL stack for propensity score
g_learner <- sl3::make_learner_stack(
  sl3::Lrnr_mean,
  sl3::Lrnr_glm_fast,
  sl3::Lrnr_earth,
  list(sl3::Lrnr_xgboost, nrounds = 50),
  list(sl3::Lrnr_xgboost, nrounds = 100),
  list(sl3::Lrnr_xgboost, nrounds = 300)
)
