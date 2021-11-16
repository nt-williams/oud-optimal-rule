library(lmtp)
library(sl3)

box::use(../R/blip, ../R/utils[...], here[here], glue[glue])

# valids options are "type1", "type2"
blip_type <- "type1"

# importing imputed data and "CATEs"
oud <- readRDS(here("data", "drv", "imputed-no-27bup.rds"))

glue("{blip_type}-blips-no27bup.rds") |>
  (\(x) here("data", "drv", x))() |>
  readRDS() -> cates

# merging blips and potential effect modifiers
.data <- purrr::map2(oud, cates, \(x, y) cbind(x[, w], y))

# creating multivariate learner stack
mv_learners <- lapply(
  list(Lrnr_xgboost$new(nrounds = 50),
       Lrnr_xgboost$new(nrounds = 100),
       Lrnr_xgboost$new(nrounds = 300),
       Lrnr_mean$new(),
       Lrnr_glm_fast$new(),
       Lrnr_earth$new()),
  \(x) sl3::make_learner(sl3::Lrnr_multivariate, x)
)

b_learner <- sl3::make_learner(sl3::Stack, mv_learners)

# estimating blips
blips <- lapply(.data, function(x) {
  if (blip_type == "type1") {
    .n <- 1:2
  }

  if (blip_type == "type2") {
    .n <- 1:3
  }

  blip$estimate_blip_multi_sl3(x, w, paste0("blip", .n), b_learner, 10)
})

# find the optimal rule
if (blip_type == "type1") {
  blips <- lapply(blips, \(x) cbind(x, 0))

  rules <- lapply(lapply(blips, \(x) x), function(x) {
    blip$find_optimal_rule(c("met", "nal", "bup"), x, TRUE) |>
      factor(levels = c("bup", "met", "nal"))
  })
}

if (blip_type == "type2") {
  rules <- lapply(lapply(blips, \(x) x), function(x) {
    blip$find_optimal_rule(c("met", "bup", "nal"), x, TRUE) |>
      factor(levels = c("bup", "met", "nal"))
  })
}

# save results
glue("estimated-super-learner-{blip_type}-blips-no27bup.rds") |>
  (\(x) here("data", "drv", x))() |>
  (\(x) saveRDS(blips, x))()

glue("estimated-super-learner-{blip_type}-rules-no27bup.rds") |>
  (\(x) here("data", "drv", x))() |>
  (\(x) saveRDS(rules, x))()
