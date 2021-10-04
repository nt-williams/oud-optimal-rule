library(lmtp)
library(purrr)

box::use(../R/rubin[...])

read_results <- \(x) readRDS(here::here("data", "drv", x))

marginal_tsms <- read_results("onestep-tsm-imputed.rds")

al_t1 <- read_results("optimal-tsms-adaptLASSO-type1.rds")
al_t2 <- read_results("optimal-tsms-adaptLASSO-type2.rds")
sl_t1 <- read_results("optimal-tsms-super-learner-type1.rds")
sl_t2 <- read_results("optimal-tsms-super-learner-type2.rds")

met <- lapply(marginal_tsms, \(x) x$met)
nal <- lapply(marginal_tsms, \(x) x$nal)
bup <- lapply(marginal_tsms, \(x) x$bup)

# TSM when receiving medication with P(A = a) = 1
rubins_rules(met)
rubins_rules(nal)
rubins_rules(bup)

# TSM under estimated optimal rules
rubins_rules(al_t1)
rubins_rules(sl_t1)
rubins_rules(al_t2)
rubins_rules(sl_t2)

# contrasts between the estimates
map(
  list(met = met, nal = nal, bup = bup),
  \(z) map2(al_t1, z, \(x, y) lmtp_contrast(x, ref = y))
) |>
  map(rubins_rules)

map(
  list(met = met, nal = nal, bup = bup),
  \(z) map2(sl_t1, z, \(x, y) lmtp_contrast(x, ref = y))
) |>
  map(rubins_rules)



