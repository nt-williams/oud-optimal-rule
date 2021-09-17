estimate_blip_binary_sl3 <- function(data, covar, blip, b_learner, nfolds) {
  folds <- origami::make_folds(data, V = nfolds)
  blip_preds <- lapply(folds, \(x) fit_blip(x, data, covar, blip, b_learner))

  purrr::map(folds, \(x) x$validation_set) |>
    purrr::reduce(c) |>
    (\(x) purrr::reduce(blip_preds, c)[order(x)])()
}

fit_blip_binary_sl3 <- function(fold, data, covar, blip, b_learner) {
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

  SL <- sl3::make_learner(
    sl3::Lrnr_sl,
    learners = b_learner,
    metalearner = sl3::make_learner("Lrnr_nnls"),
    keep_extra = FALSE
  )

  fit <- SL$train(train_task)
  fit$predict(valid_task)
}

estimate_blip_multi_sl3 <- function(data, covar, blip, b_learner, nfolds) {
  folds <- origami::make_folds(data, V = nfolds)

  blip_preds <- lapply(
    folds, \(x) fit_blip_multi_sl3(x, data, covar, blip, b_learner)
  )

  purrr::map(folds, \(x) x$validation_set) |>
    purrr::reduce(c) |>
    (\(x) purrr::reduce(blip_preds, rbind)[order(x), ])()
}

fit_blip_multi_sl3 <- function(fold, data, covar, blip, b_learner) {
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

  SL <- sl3::make_learner(
    sl3::Lrnr_sl,
    learners = b_learner,
    metalearner = sl3::make_learner(
      sl3::Lrnr_solnp,
      loss_function = sl3::loss_squared_error_multivariate,
      learner_function = sl3::metalearner_linear_multivariate
    ),
    keep_extra = FALSE
  )

  fit <- SL$train(train_task)
  fit$predict(valid_task) |> sl3::unpack_predictions()
}
