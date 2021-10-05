library(lmtp)
library(tidyverse)
library(patchwork)

box::use(../R/rubin[...], here[here])

read_results <- \(x) readRDS(here("data", "drv", x))

marginals <-
  read_results("onestep-tsm-imputed.rds") |>
  (\(r) map(c("met", "nal", "bup"), \(x) map(r, \(y) y[[x]])))()

c("optimal-tsms-adaptLASSO-type1.rds",
  "optimal-tsms-adaptLASSO-type2.rds",
  "optimal-tsms-super-learner-type1.rds",
  "optimal-tsms-super-learner-type2.rds") |>
  map(read_results) -> optimals

ans <- rbind(
  # TSM when receiving medication with P(A = a) = 1
  rubins_rules(marginals[[1]], "Methadone"),
  rubins_rules(marginals[[2]], "Naltrexone"),
  rubins_rules(marginals[[3]], "Buprenorphine"),

  # TSM under estimated optimal rules
  rubins_rules(optimals[[1]], "LASSO"),
  rubins_rules(optimals[[3]], "SL"),

  # contrasts between the estimates
  map(
    list(marginals[[1]], marginals[[2]], marginals[[3]]),
    \(z) map2(al_t1, z, \(x, y) lmtp_contrast(x, ref = y))
  ) |>
    map2_dfr(
      c("RD: methadone, LASSO",
        "RD: naltrexone, LASSO",
        "RD: buprenorphine, LASSO"),
      rubins_rules
    ),

  map(
    list(marginals[[1]], marginals[[2]], marginals[[3]]),
    \(z) map2(sl_t1, z, \(x, y) lmtp_contrast(x, ref = y))
  ) |>
    map2_dfr(
      c("RD: methadone, SL",
        "RD: naltrexone, SL",
        "RD: buprenorphine, SL"),
      rubins_rules
    )
)

ans |>
  filter(!grepl("^RD", label)) |>
  mutate(
    label = if_else(label %in% c("LASSO", "SL"), "Rule", label),
    model = c(rep("none", 3), "LASSO", "SL")
  ) |>
  ggplot(aes(
    x = factor(
      label,
      levels = c("Buprenorphine", "Methadone", "Naltrexone", "Rule")
    ),
    y = theta,
    linetype = factor(model, levels = c("none", "LASSO", "SL"))
  )) +
  geom_point(position = position_dodge(.5)) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      linetype = factor(model, levels = c("none", "LASSO", "SL"))
    ),
    width = 0.2,
    position = position_dodge(.5)
  ) +
  coord_cartesian(ylim = c(0.25, 0.75)) +
  scale_linetype_discrete(breaks = c("LASSO", "SL")) +
  labs(x = NULL,
       y = "Expected risk of relapse by 12-weeks",
       linetype = NULL) +
  theme_bw() + {
    ans |>
      filter(grepl("^RD", label)) |>
      separate(label, c("label", "model"), sep = ",") |>
      mutate(label = case_when(
        grepl("methadone", label) ~ "Methadone",
        grepl("naltrexone", label) ~ "Naltrexone",
        TRUE ~ "Buprenorphine"
      )) |>
      ggplot(aes(
        x = label,
        y = theta,
        linetype = model
      )) +
      geom_point(position = position_dodge(0.5)) +
      geom_errorbar(
        aes(
          ymin = conf.low,
          ymax = conf.high,
        ),
        width = 0.2,
        position = position_dodge(0.5)
      ) +
      scale_linetype_manual(values = c("dashed", "longdash"), guide = NULL) +
      labs(
        x = NULL,
        y = "Expected difference in risk of relapse by 12-weeks",
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(-0.35, 0)) +
      theme_bw()
  } +
  plot_layout(guides = "collect")

saveRDS(ans, here("data", "drv", "results.rds"))
