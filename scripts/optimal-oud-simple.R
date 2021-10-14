box::use(../R/blip, ../R/utils[...], here[here], glue[glue])

blip_type <- "type1"
.f <- cbind(blip1, blip2) ~ .

# importing imputed data and "CATEs"
oud <- readRDS(here("data", "drv", "imputed-coded.rds"))
cates <- readRDS(here("data", "drv", glue("{blip_type}-blips.rds")))

# merging blips and potential effect modifiers
.data <- purrr::map2(oud, cates, \(x, y) cbind(x[, w], y)) |>
  lapply(function(x) {
    model.matrix(
      ~ .,
      data = x
    )[, c("blip1", "blip2", "xrace3", "hwithdraw4", "cocdisorder", "hasEpilepsy", "ivdrug")]
  })

fit <- function(data, formula, nfolds, blip_type) {
  folds <- origami::make_folds(data, V = nfolds)

  blip_preds <- lapply(folds, function(fold) {
    .train <- origami::training(data)
    .valid <- origami::validation(data)

    fit <- lm(formula, data = as.data.frame(.train))

    list(
      coef = coef(fit),
      pred = predict(fit, as.data.frame(.valid))
    )
  })

  list(
    coef = lapply(blip_preds, \(x) x$coef),
    pred = {
      purrr::map(folds, \(x) x$validation_set) |>
        purrr::reduce(c) |>
        (\(x) purrr::reduce(lapply(blip_preds, \(x) x$pred), rbind)[order(x), ])()
    }
  )
}

# estimating blips
blips <- lapply(.data, \(x) fit(x, .f, 10, blip_type))

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
glue("estimated-lm-{blip_type}-blips.rds") |>
  (\(x) here("data", "drv", x))() |>
  (\(x) saveRDS(blips, x))()

glue("estimated-lm-{blip_type}-rules.rds") |>
  (\(x) here("data", "drv", x))() |>
  (\(x) saveRDS(rules, x))()
