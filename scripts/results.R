library(lmtp)
library(tidyverse)
library(patchwork)

box::use(../R/rubin[...], here[here])

read_results <- \(x) readRDS(here("data", "drv", x))

# marginals <-
#   read_results("onestep-tsm-imputed-no27bup.rds") |>
#   (\(r) map(c("met", "nal", "bup"), \(x) map(r, \(y) y[[x]])))()
#
# crossing("optimal-tsms-",
#          c("adaptLASSO", "super-learner"),
#          "-type1-no27bup", ".rds") |>
#   pmap_chr(paste0) |>
#   map(read_results) -> optimals
#
# ans <- rbind(
#   # TSM when receiving medication with P(A = a) = 1
#   rubins_rules(marginals[[1]], "Methadone"),
#   rubins_rules(marginals[[2]], "Naltrexone"),
#   rubins_rules(marginals[[3]], "Buprenorphine"),
#
#   # TSM under estimated optimal rules
#   rubins_rules(optimals[[1]], "LASSO"),
#   rubins_rules(optimals[[2]], "SL"),
#
#   # contrasts between the estimates
#   map(
#     list(marginals[[1]], marginals[[2]], marginals[[3]]),
#     \(z) map2(optimals[[1]], z, \(x, y) lmtp_contrast(x, ref = y))
#   ) |>
#     map2_dfr(
#       c("RD: methadone, LASSO",
#         "RD: naltrexone, LASSO",
#         "RD: buprenorphine, LASSO"),
#       rubins_rules
#     ),
#
#   map(
#     list(marginals[[1]], marginals[[2]], marginals[[3]]),
#     \(z) map2(optimals[[2]], z, \(x, y) lmtp_contrast(x, ref = y))
#   ) |>
#     map2_dfr(
#       c("RD: methadone, SL",
#         "RD: naltrexone, SL",
#         "RD: buprenorphine, SL"),
#       rubins_rules
#     ),
#
#   map(
#     list(marginals[[1]], marginals[[2]], marginals[[3]]),
#     \(z) map2(optimals[[1]], z, \(x, y) lmtp_contrast(x, ref = y, type = "rr"))
#   ) |>
#     map2_dfr(
#       c("RR: methadone, LASSO",
#         "RR: naltrexone, LASSO",
#         "RR: buprenorphine, LASSO"),
#       rubins_rules
#     ),
#
#   map(
#     list(marginals[[1]], marginals[[2]], marginals[[3]]),
#     \(z) map2(optimals[[2]], z, \(x, y) lmtp_contrast(x, ref = y, type = "rr"))
#   ) |>
#     map2_dfr(
#       c("RR: methadone, SL",
#         "RR: naltrexone, SL",
#         "RR: buprenorphine, SL"),
#       rubins_rules
#     )
# )
#
# saveRDS(ans, here::here("data", "drv", "estimates-no27bup.rds"))
ans <- readRDS(here::here("data", "drv", "estimates-no27bup.rds"))

by_study <- read_results("observed-by-study-relapse.rds") |>
  rename(label = project,
         theta = relapse,
         conf.low = "Lower bound",
         conf.high = "Upper bound")

by_medicine <- read_results("observed-by-medicine-relapse.rds") |>
  rename(label = medicine,
         theta = relapse,
         conf.low = "Lower bound",
         conf.high = "Upper bound")

ans <- full_join(ans, by_study) |>
  full_join(by_medicine) |>
  mutate(label = case_when(
    label == "bup" ~ "Obs. Buprenorphine",
    label == "met" ~ "Obs. Methadone",
    label == "nal" ~ "Obs. Naltrexone",
    TRUE ~ label
  ))

ans <- ans[c(18:20, 22, 23, 21, 1:17), ]

ragg::agg_png("plots/figure•withObs•no27bup.png", width = 8, height = 4, units = "cm", res = 400)

ans |>
  filter(!grepl("^RD|^RR", label)) |>
  mutate(
    model = case_when(
      label == "LASSO" ~ "hat(d)^lasso",
      label == "SL" ~ "hat(d)^sl",
      label %in% c("CTN0027", "CTN0030", "CTN0051", "Obs. Buprenorphine", "Obs. Methadone", "Obs. Naltrexone") ~ "Observed",
      TRUE ~ "TSM"
    ),
    label = case_when(
      label %in% c("LASSO", "SL") ~ "Rule",
      label == "Obs. Buprenorphine" ~ "Buprenorphine",
      label == "Obs. Methadone" ~ "Methadone",
      label == "Obs. Naltrexone" ~ "Naltrexone",
      TRUE ~ label
    )
  ) |>
  ggplot(aes(
    x = factor(
      label,
      levels = c("CTN0027", "CTN0030", "CTN0051", "Buprenorphine", "Methadone", "Naltrexone", "Rule")
    ),
    y = theta,
    linetype = factor(model, levels = c("Observed", "TSM", "hat(d)^lasso", "hat(d)^sl"))
  )) +
  geom_point(position = position_dodge(.75), size = 0.15) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      linetype = factor(model, levels = c("Observed", "TSM", "hat(d)^lasso", "hat(d)^sl"))
    ),
    width = 0.15,
    position = position_dodge(.75),
    size = 0.15
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  coord_cartesian(ylim = c(0.1, 0.7)) +
  scale_linetype_manual(
    breaks = c("Observed", "TSM", "hat(d)^lasso", "hat(d)^sl"),
    labels = scales::parse_format(),
    values = c(1, 2, 3, 4)
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
    ans |>
      filter(grepl("^RD", label)) |>
      separate(label, c("label", "model"), sep = ",") |>
      mutate(
        label = case_when(
          grepl("methadone", label) ~ "Methadone",
          grepl("naltrexone", label) ~ "Naltrexone",
          TRUE ~ "Buprenorphine"
        ),
        model = if_else(model == " LASSO", "hat(d)^lasso", "hat(d)^sl")
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
        width = 0.2,
        position = position_dodge(0.75),
        size = 0.15
      ) +
      geom_hline(yintercept = 0, size = 0.15) +
      scale_linetype_manual(values = c(3, 4), guide = NULL) +
      labs(
        x = NULL,
        y = "Expected difference in risk of relapse by 12-weeks",
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(-0.5, 0.1)) +
      theme_light(base_size = 3,
               base_line_size = 0.2,
               base_rect_size = 0.2) +
      theme(legend.text = element_text(size = 3),
            legend.key.size = unit(0.2, "cm"),
            axis.ticks.x = element_blank()) +
      guides(shape = guide_legend(override.aes = list(size = 0.5)))
  } +
  plot_layout(guides = "collect", widths = c(0.65, 0.35))

dev.off()
