# Currently only considering main effects...

library(lmtp)
library(glmnet)

box::use(../R/blip, ../R/utils[...], here[here], glue[glue])

blip_type <- "type1"
crossfit <- FALSE

V <- ifelse(crossfit, 10, 1)

# importing imputed data and "CATEs"
oud <- readRDS(here("data", "drv", "imputed-no-27bup.rds"))
cates <- readRDS(here("data", "drv", glue("{blip_type}-blips-no27bup.rds")))

# merging blips and potential effect modifiers
.data <- purrr::map2(oud, cates, \(x, y) cbind(x[, w], y))

min_weight <- function(x) {
  .min <- max.col(x * -1)
  out <- vector("double", nrow(x))
  for (i in 1:nrow(x)) {
    out[i] <- x[i, .min[i]]
  }
  out
}

adaptive_lasso <- function(data, covar, nfolds, blip_type) {
  folds <- origami::make_folds(data, V = nfolds)

  if (nfolds == 1) {
    folds[[1]]$training_set <- folds[[1]]$validation_set
  }

  blip_preds <- lapply(folds, \(x) fit_adaptive_lasso(x, data, covar, blip_type))

  list(
    fit = lapply(blip_preds, \(x) x$fit),
    selected = lapply(blip_preds, \(x) x$selected),
    pred = {
      purrr::map(folds, \(x) x$validation_set) |>
        purrr::reduce(c) |>
        (\(x) purrr::reduce(lapply(blip_preds, \(x) x$pred), rbind)[order(x), ])()
    }
  )
}

fit_adaptive_lasso <- function(fold, data, covar, type) {
  .train <- origami::training(data)
  .valid <- origami::validation(data)

  if (type == "type1") {
    .f <- as.formula("cbind(blip1, blip2) ~ .")
    .n <- 1:2
  }

  if (type == "type2") {
    .f <- as.formula("cbind(blip1, blip2, blip3) ~ .")
    .n <- 1:3
  }

  penalty <-
    lm(.f, data = .train) |>
    (\(x) 1 / abs(x$coefficients))() |>
    min_weight() |>
    (\(x) ifelse(is.na(x), 2000, x))()

  fit <- cv.glmnet(
    model.matrix(~ ., .train[, covar])[, -1],
    as.matrix(.train[, paste0("blip", .n)]),
    family = "mgaussian", penalty.factor = penalty[-1]
  )

  list(
    fit = fit,
    selected = coef(fit, s = "lambda.min"),
    pred = predict(
      fit, s = "lambda.min", gamma = c(1), relax = TRUE,
      newx = model.matrix(~ ., .valid[, covar])[, -1]
    )[, , 1]
  )
}

# estimating blips
blips <- lapply(.data, \(x) adaptive_lasso(x, w, V, blip_type))

# find the optimal rule
if (blip_type == "type1") {
  blips <- lapply(blips, function(x) {
    x$pred <- cbind(x$pred, 0)
    x
  })

  rules <- lapply(lapply(blips, \(x) x$pred), function(x) {
    blip$find_optimal_rule(c("met", "nal", "bup"), x, TRUE) |>
      factor(levels = c("bup", "met", "nal"))
  })
}

if (blip_type == "type2") {
  rules <- lapply(lapply(blips, \(x) x$pred), function(x) {
    blip$find_optimal_rule(c("met", "bup", "nal"), x, TRUE) |>
      factor(levels = c("bup", "met", "nal"))
  })
}

# save results
glue("estimated-adaptLASSO-{crossfit}-{blip_type}-blips-no27bup.rds") |>
  (\(x) here("data", "drv", x))() |>
  (\(x) saveRDS(blips, x))()

glue("estimated-adaptLASSO-{crossfit}-{blip_type}-rules-no27bup.rds") |>
  (\(x) here("data", "drv", x))() |>
  (\(x) saveRDS(rules, x))()
