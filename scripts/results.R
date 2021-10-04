library(lmtp)
library(purrr)

box::use(../R/rubin[...], here[here])

read_results <- \(x) readRDS(here("data", "drv", x))

marginal_tsms <- read_results("onestep-tsm-imputed.rds")

al_t1 <- read_results("optimal-tsms-adaptLASSO-type1.rds")
al_t2 <- read_results("optimal-tsms-adaptLASSO-type2.rds")
sl_t1 <- read_results("optimal-tsms-super-learner-type1.rds")
sl_t2 <- read_results("optimal-tsms-super-learner-type2.rds")

met <- lapply(marginal_tsms, \(x) x$met)
nal <- lapply(marginal_tsms, \(x) x$nal)
bup <- lapply(marginal_tsms, \(x) x$bup)

ans <- list(
  # TSM when receiving medication with P(A = a) = 1
  all_methadone = rubins_rules(met),
  all_naltrexone = rubins_rules(nal),
  all_bupenorphine = rubins_rules(bup),

  # TSM under estimated optimal rules
  opt_lasso_t1 = rubins_rules(al_t1),
  opt_sl_t1 = rubins_rules(sl_t1),
  opt_lasso_t2 = rubins_rules(al_t2),
  opt_sl_t2 = rubins_rules(sl_t2),

  # contrasts between the estimates
  contrast_lasso_t1 = map(
    list(met = met, nal = nal, bup = bup),
    \(z) map2(al_t1, z, \(x, y) lmtp_contrast(x, ref = y))
  ) |>
    map(rubins_rules),

  contrast_sl_t1 = map(
    list(met = met, nal = nal, bup = bup),
    \(z) map2(sl_t1, z, \(x, y) lmtp_contrast(x, ref = y))
  ) |>
    map(rubins_rules)
)

saveRDS(ans, here("data", "drv", "results.rds"))
