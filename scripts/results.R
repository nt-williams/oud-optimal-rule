library(lmtp)
library(tidyverse)
library(patchwork)

box::use(../R/rubin[...], here[here])

read_results <- \(x) readRDS(here("data", "drv", x))

marginals <-
  read_results("onestep-tsm-imputed-no27bup.rds") |>
  (\(r) map(c("met", "nal", "bup"), \(x) map(r, \(y) y[[x]])))()

crossing("optimal-tsms-",
         c("adaptLASSO", "super-learner"),
         "-type1-no27bup", ".rds") |>
  pmap_chr(paste0) |>
  map(read_results) -> optimals

ans <- rbind(
  # TSM when receiving medication with P(A = a) = 1
  rubins_rules(marginals[[1]], "Methadone"),
  rubins_rules(marginals[[2]], "Naltrexone"),
  rubins_rules(marginals[[3]], "Buprenorphine"),

  # TSM under estimated optimal rules
  rubins_rules(optimals[[1]], "LASSO"),
  rubins_rules(optimals[[2]], "SL"),

  # contrasts between the estimates
  map(
    list(marginals[[1]], marginals[[2]], marginals[[3]]),
    \(z) map2(optimals[[1]], z, \(x, y) lmtp_contrast(x, ref = y))
  ) |>
    map2_dfr(
      c("RD: methadone, LASSO",
        "RD: naltrexone, LASSO",
        "RD: buprenorphine, LASSO"),
      rubins_rules
    ),

  map(
    list(marginals[[1]], marginals[[2]], marginals[[3]]),
    \(z) map2(optimals[[2]], z, \(x, y) lmtp_contrast(x, ref = y))
  ) |>
    map2_dfr(
      c("RD: methadone, SL",
        "RD: naltrexone, SL",
        "RD: buprenorphine, SL"),
      rubins_rules
    ),

  map(
    list(marginals[[1]], marginals[[2]], marginals[[3]]),
    \(z) map2(optimals[[1]], z, \(x, y) lmtp_contrast(x, ref = y, type = "rr"))
  ) |>
    map2_dfr(
      c("RR: methadone, LASSO",
        "RR: naltrexone, LASSO",
        "RR: buprenorphine, LASSO"),
      rubins_rules
    ),

  map(
    list(marginals[[1]], marginals[[2]], marginals[[3]]),
    \(z) map2(optimals[[2]], z, \(x, y) lmtp_contrast(x, ref = y, type = "rr"))
  ) |>
    map2_dfr(
      c("RR: methadone, SL",
        "RR: naltrexone, SL",
        "RR: buprenorphine, SL"),
      rubins_rules
    )
)

saveRDS(ans, here::here("data", "drv", "estimates-no27bup.rds"))

ans |>
  filter(!grepl("^RD|^RR", label)) |>
  mutate(
    label = if_else(label %in% c("LASSO", "SL"), "Rule", label),
    model = c(rep("none", 3), "hat(d)(v)^lasso", "hat(d)(v)^sl")
  ) |>
  ggplot(aes(
    x = factor(
      label,
      levels = c("Buprenorphine", "Methadone", "Naltrexone", "Rule")
    ),
    y = theta,
    linetype = factor(model, levels = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl"))
  )) +
  geom_point(position = position_dodge(.75)) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      linetype = factor(model, levels = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl"))
    ),
    width = 0.2,
    position = position_dodge(.75)
  ) +
  coord_cartesian(ylim = c(0.2, 0.75)) +
  scale_linetype_discrete(
    breaks = c("hat(d)(v)^lasso", "hat(d)(v)^sl"),
    labels = scales::parse_format()
  ) +
  labs(
    x = NULL,
    y = "Expected risk of relapse by 12-weeks",
    linetype = NULL
  ) +
  theme_bw() +
  theme(legend.text.align = 0) + {
    ans |>
      filter(grepl("^RD", label)) |>
      separate(label, c("label", "model"), sep = ",") |>
      mutate(
        label = case_when(
          grepl("methadone", label) ~ "Methadone",
          grepl("naltrexone", label) ~ "Naltrexone",
          TRUE ~ "Buprenorphine"
        ),
        model = if_else(model == " LASSO", "hat(d)(v)^lasso", "hat(d)(v)^sl")
      ) |>
      ggplot(aes(
        x = label,
        y = theta,
        linetype = model
      )) +
      geom_point(position = position_dodge(0.75)) +
      geom_errorbar(
        aes(
          ymin = conf.low,
          ymax = conf.high,
        ),
        width = 0.2,
        position = position_dodge(0.75)
      ) +
      geom_hline(yintercept = 0) +
      scale_linetype_manual(values = c("22", "42"), guide = NULL) +
      labs(
        x = NULL,
        y = "Expected difference in risk of relapse by 12-weeks",
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(-0.3, 0.1)) +
      theme_bw()
  } +
  plot_layout(guides = "collect")

ggsave(here("plots", "figure-no27bup.png"))
