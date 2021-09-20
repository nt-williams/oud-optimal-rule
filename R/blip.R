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

estimate_blip_multi_singlevar <- function(data, covar, blip, nfolds) {
  folds <- origami::make_folds(data, V = nfolds)

  blip_preds <- lapply(
    folds, \(x) fit_blip_multi_singlevar(x, data, covar, blip)
  )

  list(
    preds = {
      purrr::map(folds, \(x) x$validation_set) |>
        purrr::reduce(c) |>
        (\(x) purrr::reduce(lapply(blip_preds, \(x) x$preds), rbind)[order(x), ])()
    },
    selected = lapply(blip_preds, \(x) x$selected)
  )
}

fit_blip_multi_singlevar <- function(fold, data, covar, blip) {
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

  corRank_screen_then_glm <- sl3::make_learner(
    sl3::Lrnr_multivariate,
    Pipeline$new(
      Lrnr_pkg_SuperLearner_screener$new("screen.corRankBest"),
      Lrnr_glm_fast$new()
    )
  )

  fit <- corRank_screen_then_glm$train(train_task)

  list(
    preds = fit$predict(valid_task) |> sl3::unpack_predictions(),
    selected = lapply(fit$fit_object$outcome_fits, get_selected)
  )
}

get_selected <- function(fit) {
  fit$fit_object$learner_fits$Lrnr_pkg_SuperLearner_screener_screen.corRankBest$fit_object$selected
}

screen.corRankBest <- function(Y, X, family, method = "kendall", rank = 1, ...) {
  listp <- apply(X, 2, function(x, Y, method) {
    ifelse(stats::var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
  }, Y = Y, method = method)

  rank(listp) <= rank
}

estimate_blip_multi_lasso <- function(data, covar, blip, nfolds) {

}

fit_blip_multi_lasso <- function(fold, data, covar, blip) {

}

find_optimal_rule <- function(levels, blips, minimize = FALSE) {
  if (minimize) {
    blips <- blips * -1
  }
  levels[max.col(blips)]
}

type_2_blip <- function(...) {
  objs <- list(...)
  purrr::map_dfc(objs, \(x) x$eif) |>
    stats::setNames(paste0("blip", seq_along(objs))) |>
    as.data.frame() |>
    (\(x) x - rowMeans(x))()
}
