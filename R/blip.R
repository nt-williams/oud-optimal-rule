estimate_blip <- function(data, covar, blip, b_learner, nfolds) {
  folds <- origami::make_folds(data, V = nfolds)
  blip_preds <- lapply(folds, \(x) fit_blip(x, data, covar, blip, b_learner))

  purrr::map(folds, \(x) x$validation_set) |>
    purrr::reduce(c) |>
    (\(x) purrr::reduce(blip_preds, c)[order(x)])()
}

fit_blip <- function(fold, data, covar, blip, b_learner) {
  train <- origami::training(data, fold)
  valid <- origami::validation(data, fold)

  train_task <- sl3::sl3_Task$new(
    data = train,
    covariates = covar,
    outcome = blip,
    outcome_type = "continuous"
  )

  valid_task <- sl3::sl3_Task$new(
    data = valid,
    covariates = covar,
    outcome_type = "continuous"
  )

  stack <- sl3::make_learner(
    Lrnr_sl,
    learners = b_learner,
    metalearner = sl3::make_learner("Lrnr_nnls"),
    keep_extra = FALSE
  )

  fit <- stack$train(train_task)
  fit$predict(valid_task)
}
