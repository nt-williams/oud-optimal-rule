box::use(../R/blip, ../R/utils[...])

blips <- readRDS(here::here("data", "drv", "estimated-adaptLASSO-type2-blips.rds"))

rules <- lapply(lapply(blips, \(x) x$pred), function(x) {
  blip$find_optimal_rule(c("met", "bup", "nal"), x, TRUE) |>
    factor(levels = c("bup", "met", "nal"))
})

saveRDS(rules, here::here("data", "drv", "estimated-adaptLASSO-type2-rules.rds"))
