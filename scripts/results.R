library(lmtp)

box::use(../R/rubin[...])

read_results <- \(x) readRDS(here::here("data", "drv", x))

marginal_tsms <- read_results("onestep-tsm-imputed.rds")

al_t1 <- read_results("optimal-tsms-adaptLASSO-type1.rds")
al_t2 <- read_results("optimal-tsms-adaptLASSO-type2.rds")
sl_t1 <- read_results("optimal-tsms-super-learner-type1.rds")
sl_t2 <- read_results("optimal-tsms-super-learner-type2.rds")

rubins_rules(lapply(marginal_tsms, \(x) x$met))
rubins_rules(lapply(marginal_tsms, \(x) x$nal))
rubins_rules(lapply(marginal_tsms, \(x) x$bup))
rubins_rules(al_t1)
rubins_rules(al_t2)
rubins_rules(sl_t1)
rubins_rules(sl_t2)
