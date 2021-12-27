library(lmtp, lib = "~/Desktop/lmtp-alternates/sl3")
library(tidyverse)
library(patchwork)
library(here)

source("R/rubin.R")

oud <- readRDS(here::here("data", "drv", "imputed-no-27bup.rds"))

oud |>
  (\(x) x[[1]])() |>
  (\(x) x$project == "27")() |>
  (\(x) oud[[1]][x, "week_12_relapse"])() -> obs_methadone

se = sqrt(var(obs_methadone - mean(obs_methadone)) / length(obs_methadone))

obs <- data.frame(
  label = "Obs. Methadone",
  theta = mean(obs_methadone),
  se = se,
  alpha = 0.05,
  conf.low = mean(obs_methadone) - qnorm(0.975) * se,
  conf.high = mean(obs_methadone) + qnorm(0.975) * se
)

read_results <- \(x) readRDS(here("data", "drv", x))

marginals <-
  read_results("onestep-tsm-imputed-no27.rds") |>
  (\(r) map(c("nal", "bup"), \(x) map(r, \(y) y[[x]])))()

names(marginals) <- c("nal", "bup")

crossing("optimal-tsms-",
         c("adaptLASSO", "super-learner"),
         "-type1-no27bup-projectSubsets", ".rds") |>
  pmap_chr(paste0) |>
  map(read_results) -> optimals

names(optimals) <- c("lasso", "sl")

make_lmtp <- function(eif) {
  out <- list(
    estimator = "SDR",
    theta = mean(eif),
    id = 1:length(eif),
    eif = eif,
    outcome_type = "binomial"
  )
  class(out) <- "lmtp"
  out
}

in27 <- rbind(
  obs,

  # TSM under estimated optimal rules
  rubins_rules(optimals[["lasso"]][["in27_eif"]], "LASSO"),
  rubins_rules(optimals[["sl"]][["in27_eif"]], "SL"),

  # contrasts between the estimates
  map(optimals[["lasso"]][["in27_eif"]], \(x) lmtp_contrast(make_lmtp(x), ref = make_lmtp(obs_methadone))) |>
    rubins_rules("RD: methadone, LASSO"),

  map(optimals[["sl"]][["in27_eif"]], \(x) lmtp_contrast(make_lmtp(x), ref = make_lmtp(obs_methadone))) |>
    rubins_rules("RD: methadone, SL"),

  map(optimals[["lasso"]][["in27_eif"]], \(x) lmtp_contrast(make_lmtp(x), ref = make_lmtp(obs_methadone), type = "rr")) |>
    rubins_rules("RR: methadone, LASSO"),

  map(optimals[["sl"]][["in27_eif"]], \(x) lmtp_contrast(make_lmtp(x), ref = make_lmtp(obs_methadone), type = "rr")) |>
    rubins_rules("RR: methadone, SL")
)

no27 <- rbind(
  # TSM when receiving medication with P(A = a) = 1
  rubins_rules(marginals[["nal"]], "Naltrexone"),
  rubins_rules(marginals[["bup"]], "Buprenorphine"),

  # TSM under estimated optimal rules
  rubins_rules(optimals[["lasso"]][["no27_eif"]], "LASSO"),
  rubins_rules(optimals[["sl"]][["no27_eif"]], "SL"),

  # contrasts between the estimates
  map(
    list(marginals[["nal"]], marginals[["bup"]]),
    \(z) map2(optimals[["lasso"]][["no27_eif"]], z, \(x, y) lmtp_contrast(make_lmtp(x), ref = y))
  ) |>
    map2_dfr(
      c("RD: naltrexone, LASSO",
        "RD: buprenorphine, LASSO"),
      rubins_rules
    ),

  map(
    list(marginals[["nal"]], marginals[["bup"]]),
    \(z) map2(optimals[["sl"]][["no27_eif"]], z, \(x, y) lmtp_contrast(make_lmtp(x), ref = y))
  ) |>
    map2_dfr(
      c("RD: naltrexone, SL",
        "RD: buprenorphine, SL"),
      rubins_rules
    ),

  map(
    list(marginals[["nal"]], marginals[["bup"]]),
    \(z) map2(optimals[["lasso"]][["no27_eif"]], z, \(x, y) lmtp_contrast(make_lmtp(x), ref = y, type = "rr"))
  ) |>
    map2_dfr(
      c("RR: naltrexone, LASSO",
        "RR: buprenorphine, LASSO"),
      rubins_rules
    ),

  map(
    list(marginals[["nal"]], marginals[["bup"]]),
    \(z) map2(optimals[["sl"]][["no27_eif"]], z, \(x, y) lmtp_contrast(make_lmtp(x), ref = y, type = "rr"))
  ) |>
    map2_dfr(
      c("RR: naltrexone, SL",
        "RR: buprenorphine, SL"),
      rubins_rules
    )
)

ragg::agg_png("plots/figure•no27.png", width = 8, height = 3.5, units = "cm", res = 400)

no27 |>
  filter(!grepl("^RD|^RR", label)) |>
  mutate(
    label = if_else(label %in% c("LASSO", "SL"), "Rule", label),
    model = c(rep("none", 2), "hat(d)(v)^lasso", "hat(d)(v)^sl")
  ) |>
  ggplot(aes(
    x = factor(
      label,
      levels = c("Buprenorphine", "Naltrexone", "Rule")
    ),
    y = theta,
    linetype = factor(model, levels = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl"))
  )) +
  geom_point(position = position_dodge(.75), size = 0.15) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      linetype = factor(model, levels = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl"))
    ),
    width = 0.15,
    size = 0.15,
    position = position_dodge(.75)
  ) +
  coord_cartesian(ylim = c(0.1, 0.7)) +
  scale_linetype_discrete(
    breaks = c("hat(d)(v)^lasso", "hat(d)(v)^sl"),
    labels = scales::parse_format()
  ) +
  labs(
    x = NULL,
    y = "Expected risk of relapse by 12-weeks",
    linetype = NULL
  ) +
  theme_light(base_size = 3,
              base_line_size = 0.2,
              base_rect_size = 0.2) +
  theme(legend.text.align = 0,
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2, "cm"),
        axis.ticks.x = element_blank()) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) + {
    no27 |>
      filter(grepl("^RD", label)) |>
      separate(label, c("label", "model"), sep = ",") |>
      mutate(
        label = case_when(
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
      geom_point(position = position_dodge(0.75), size = 0.15) +
      geom_errorbar(
        aes(
          ymin = conf.low,
          ymax = conf.high,
        ),
        width = 0.15,
        position = position_dodge(0.75),
        size = 0.15
      ) +
      geom_hline(yintercept = 0, size = 0.15) +
      scale_linetype_manual(values = c("22", "42"), guide = NULL) +
      labs(
        x = NULL,
        y = "Expected difference in risk of relapse by 12-weeks",
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(-0.5, 0.1)) +
      theme_light(base_size = 3,
                  base_line_size = 0.2,
                  base_rect_size = 0.2) +
      theme(legend.text.align = 0,
            legend.text = element_text(size = 3),
            legend.key.size = unit(0.2, "cm"),
            axis.ticks.x = element_blank()) +
      guides(shape = guide_legend(override.aes = list(size = 0.5)))
  } +
  plot_layout(guides = "collect")

dev.off()

ragg::agg_png("plots/figure•in27.png", width = 8, height = 3.5, units = "cm", res = 400)

in27 |>
  filter(!grepl("^RD|^RR", label)) |>
  mutate(
    label = if_else(label %in% c("LASSO", "SL"), "Rule", label),
    model = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl")
  ) |>
  ggplot(aes(
    x = factor(
      label,
      levels = c("Obs. Methadone", "Rule")
    ),
    y = theta,
    linetype = factor(model, levels = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl"))
  )) +
  geom_point(position = position_dodge(.75), size = 0.15) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      linetype = factor(model, levels = c("none", "hat(d)(v)^lasso", "hat(d)(v)^sl"))
    ),
    width = 0.15,
    position = position_dodge(.75),
    size = 0.15
  ) +
  coord_cartesian(ylim = c(0.1, 0.7)) +
  scale_linetype_discrete(
    breaks = c("hat(d)(v)^lasso", "hat(d)(v)^sl"),
    labels = scales::parse_format()
  ) +
  labs(
    x = NULL,
    y = "Expected risk of relapse by 12-weeks",
    linetype = NULL
  ) +
  theme_light(base_size = 3,
              base_line_size = 0.2,
              base_rect_size = 0.2) +
  theme(legend.text.align = 0,
        legend.text = element_text(size = 3),
        legend.key.size = unit(0.2, "cm"),
        axis.ticks.x = element_blank()) +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) + {
    in27 |>
      filter(grepl("^RD", label)) |>
      separate(label, c("label", "model"), sep = ",") |>
      mutate(
        label = "Obs. Methadone",
        model = if_else(model == " LASSO", "hat(d)(v)^lasso", "hat(d)(v)^sl")
      ) |>
      ggplot(aes(
        x = label,
        y = theta,
        linetype = model
      )) +
      geom_point(position = position_dodge(0.75), size = 0.15) +
      geom_errorbar(
        aes(
          ymin = conf.low,
          ymax = conf.high,
        ),
        width = 0.15,
        position = position_dodge(0.75),
        size = 0.15
      ) +
      geom_hline(yintercept = 0, size = 0.15) +
      scale_linetype_manual(values = c("22", "42"), guide = NULL) +
      labs(
        x = NULL,
        y = "Expected difference in risk of relapse by 12-weeks",
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(-0.5, 0.1)) +
      theme_light(base_size = 3,
                  base_line_size = 0.2,
                  base_rect_size = 0.2) +
      theme(legend.text.align = 0,
            legend.text = element_text(size = 3),
            legend.key.size = unit(0.2, "cm"),
            axis.ticks.x = element_blank()) +
      guides(shape = guide_legend(override.aes = list(size = 0.5)))
  } +
  plot_layout(guides = "collect")

dev.off()
