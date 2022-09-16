library(lmtp)
library(tidyverse)
library(patchwork)

box::use(../R/rubin[...], here[here])

read_results <- \(x) readRDS(here("data", "drv", x))
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

pt_to_mm <- function(pt) {
  pt / ggplot2::.pt
}

p1 <- ans |>
  filter(!grepl("^RD|^RR", label)) |>
  mutate(
    model = case_when(
      label == "LASSO" ~ "hat(d)^lasso",
      label == "SL" ~ "hat(d)^sl",
      label %in% c("CTN0027", "CTN0030", "CTN0051", "Obs. Buprenorphine", "Obs. Methadone", "Obs. Naltrexone") ~ "Observed",
      TRUE ~ "Estimated"
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
    linetype = factor(model, levels = c("Observed", "Estimated", "hat(d)^lasso", "hat(d)^sl"))
  )) +
  geom_point(position = position_dodge(.75), size = pt_to_mm(1)) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      linetype = factor(model, levels = c("Observed", "Estimated", "hat(d)^lasso", "hat(d)^sl"))
    ),
    width = 0.15,
    position = position_dodge(.75),
    size = pt_to_mm(1)
  ) +
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  coord_cartesian(ylim = c(0.2, 1)) +
  scale_linetype_manual(
    breaks = c("Observed", "Estimated", "hat(d)^lasso", "hat(d)^sl"),
    labels = scales::parse_format(),
    values = c(1, 2, 3, 4)
  ) +
  labs(
    x = "Treatment",
    y = "Expected Risk of Relapse by 12-weeks",
    linetype = expression(underline("Estimand"))
  ) +
  theme_classic(base_size = 8,
                base_line_size = 0.4,
                base_rect_size = 0.4) +
  theme(legend.position = c(0.15, 0.775),
        legend.background = element_rect(fill = "white",
                                         color = "black",
                                         size = pt_to_mm(1)),
        text = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.line = element_line(size = pt_to_mm(1)),
        axis.ticks = element_line(size = pt_to_mm(1)),
        legend.text = element_text(size = 9),
        legend.text.align = 0,
        legend.title.align = 0.5) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)))

p2 <- ans |>
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
  geom_point(position = position_dodge(0.75), size = pt_to_mm(1)) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high,
    ),
    width = 0.2,
    position = position_dodge(0.75),
    size = pt_to_mm(1)
  ) +
  geom_hline(yintercept = 0, size = pt_to_mm(1), linetype = "dashed") +
  scale_linetype_manual(values = c(3, 4), guide = NULL) +
  scale_x_discrete(guide = guide_axis(angle = -45)) +
  scale_y_continuous(breaks = c(0, -0.2, -0.4),
                     labels = c("0", "-0.2", "-0.4")) +
  labs(
    x = "Treatment",
    y = "Expected Difference in Risk \nof Relapse by 12-weeks",
    linetype = NULL
  ) +
  coord_cartesian(ylim = c(-0.5, 0.1)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 9),
        axis.text = element_text(size = 9, colour = "black"),
        axis.line = element_line(size = pt_to_mm(1)),
        axis.ticks = element_line(size = pt_to_mm(1)))

pdf.options(encoding='ISOLatin2.enc')

pdf("plots/AJE-00224-2022 Rudolph Figure 1.pdf", width = 8, height = 4)
p1 + p2 + plot_layout(widths = c(0.65, 0.35))
dev.off()

ggsave("plots/AJE-00224-2022 Rudolph Figure 1A.eps",
       device = "eps",
       plot = p1,
       width = .65*8, height = 4, units = "in")

ggsave("plots/AJE-00224-2022 Rudolph Figure 1B.eps",
       device = "eps",
       plot = p2,
       width = .35*8, height = 4, units = "in")
