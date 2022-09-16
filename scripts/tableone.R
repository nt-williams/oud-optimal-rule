library(tidyverse)

box::use(../R/utils[...])

local({
  load(here::here("data", "src", "clean_combined_imputed_data9-8-21.Rdata"))
  baseline <<- ALT_patients_with_outcomes_02
})

oud <- readRDS(here::here("data", "drv", "imputed-no-27bup.rds"))[[1]]

baseline <-
  select(baseline, who, project, all_of(w)) |>
  mutate(
    across(starts_with("has"), ~ as.numeric(.x == "yes")),
    across(ends_with("disorder"), ~ as.numeric(.x == "yes")),
    across(c("bamphetamine30_base", "bcannabis30_base", "bbenzo30_base", "ivdrug"),
           ~ as.numeric(.x == "yes")),
    who = as.character(who)
  )

dat <-
  select(oud, who, weeks_to_relapse, week_12_relapse) |>
  left_join(baseline) |>
  as_tibble()

pmean <- \(x, ...) sprintf("%0.1f\\%%", mean(x, na.rm = TRUE) * 100)
cmean <- \(x, ...) sprintf("%.2f", mean(x, na.rm = TRUE))
csd <- \(x, ...) sprintf("%.2f", sd(x, na.rm = TRUE))

stats <- function(data, ...) {
  summarise(data,
            n = n(),
            women = pmean(sex == "female"),
            age_sd = csd(age),
            age = cmean(age),
            race_1 = pmean(xrace == "1"),
            race_2 = pmean(xrace == "2"),
            race_3 = pmean(xrace == "3"),
            race_4 = pmean(xrace == "4"),
            iv = pmean(ivdrug),
            alc = pmean(alcdisorder),
            coc = pmean(cocdisorder),
            brain_damage = pmean(hasBrainDamage),
            epilepsy = pmean(hasEpilepsy),
            schiz = pmean(hasSchiz),
            bipolar = pmean(hasBipolar),
            anxiety = pmean(hasAnxPan),
            cannabis = pmean(bcannabis30_base),
            meth = pmean(bamphetamine30_base),
            benzo = pmean(bbenzo30_base),
            hwithdraw_sd = csd(as.numeric(hwithdraw)),
            hwithdraw = cmean(as.numeric(hwithdraw)),
            rweek_sd = csd(weeks_to_relapse),
            rweek = cmean(weeks_to_relapse),
            w12 = sum(week_12_relapse)
  )
}

total <- stats(dat)

by_project <-
  group_by(dat, project) |>
  nest() |>
  mutate(stats = map(data, \(x) stats(x))) |>
  (\(x) x$stats)() |>
  setNames(c("27", "51", "30"))

# Produces LaTeX for table 1, written to the clipboard
glue::glue(
  "\\begin{table}
  \\caption{Clinical characteristics of patients by CTN trial.}
  \\centering\\footnotesize
  \\begin{tabular}[t]{lcccc}
  \\toprule
  & All & CTN0027 & CTN0030 & CTN0051 \\\\ \
  \\midrule
  \\addlinespace[0.3em]
  \\midrule
  N & <total$n> & <by_project$`27`$n> & <by_project$`30`$n> & <by_project$`51`$n> \\\\ \
  Age & <total$age> (<total$age_sd>) & <by_project$`27`$age> (<by_project$`27`$age_sd>) & <by_project$`30`$age> (<by_project$`30`$age_sd>) & <by_project$`51`$age> (<by_project$`51`$age_sd>) \\\\ \
  Women & <total$women> & <by_project$`27`$women> & <by_project$`30`$women> & <by_project$`51`$women> \\\\ \
  \\multicolumn{6}{l}{Race/ethnicity} \\\\ \
  \\hspace{1em} Non-Hispanic white & <total$race_1> & <by_project$`27`$race_1> & <by_project$`30`$race_1> & <by_project$`51`$race_1> \\\\ \
  \\hspace{1em} Non-Hispanic Black & <total$race_2> & <by_project$`27`$race_2> & <by_project$`30`$race_2> & <by_project$`51`$race_2> \\\\ \
  \\hspace{1em} Hispanic & <total$race_3> & <by_project$`27`$race_3> & <by_project$`30`$race_3> & <by_project$`51`$race_3> \\\\ \
  \\hspace{1em} Other (including multiracial) & <total$race_4> & <by_project$`27`$race_4> & <by_project$`30`$race_4> & <by_project$`51`$race_4> \\\\ \
  Current IV drug use & <total$iv> & <by_project$`27`$iv> & <by_project$`30`$iv> & <by_project$`51`$iv> \\\\ \
  Current cannabis use & <total$cannabis> & <by_project$`27`$cannabis> & <by_project$`30`$cannabis> & <by_project$`51`$cannabis> \\\\ \
  Current amphetamine use & <total$meth> & <by_project$`27`$meth> & <by_project$`30`$meth> & <by_project$`51`$meth> \\\\ \
  Current benzodiazepine drug use & <total$benzo> & <by_project$`27`$benzo> & <by_project$`30`$benzo> & <by_project$`51`$benzo> \\\\ \
  Alcohol use disorder & <total$alc> & <by_project$`27`$alc> & <by_project$`30`$alc> & <by_project$`51`$alc> \\\\ \
  Cocaine use disorder & <total$coc> & <by_project$`27`$coc> & <by_project$`30`$coc> & <by_project$`51`$coc> \\\\ \
  Neurological injury & <total$brain_damage> & <by_project$`27`$brain_damage> & <by_project$`30`$brain_damage> & <by_project$`51`$brain_damage> \\\\ \
  History of epilepsy & <total$epilepsy> & <by_project$`27`$epilepsy> & <by_project$`30`$epilepsy> & <by_project$`51`$epilepsy> \\\\ \
  History of schizophrenia & <total$schiz> & <by_project$`27`$schiz> & <by_project$`30`$schiz> & <by_project$`51`$schiz> \\\\ \
  History of bipolar disorder & <total$bipolar> & <by_project$`27`$bipolar> & <by_project$`30`$bipolar> & <by_project$`51`$bipolar> \\\\ \
  History of anxiety disorder & <total$anxiety> & <by_project$`27`$anxiety> & <by_project$`30`$anxiety> & <by_project$`51`$anxiety> \\\\ \
  Opioid withdrawl discomfort (1-4) & <total$hwithdraw> (<total$hwithdraw_sd>) & <by_project$`27`$hwithdraw> (<by_project$`27`$hwithdraw_sd>) & <by_project$`30`$hwithdraw> (<by_project$`30`$hwithdraw_sd>) & <by_project$`51`$hwithdraw> (<by_project$`51`$hwithdraw_sd>) \\\\ \
  Week of relapse & <total$rweek> (<total$rweek_sd>) & <by_project$`27`$rweek> (<by_project$`27`$rweek_sd>) & <by_project$`30`$rweek> (<by_project$`30`$rweek_sd>) & <by_project$`51`$rweek> (<by_project$`51`$rweek_sd>) \\\\ \
  No. of relapses by week 12 & <total$w12> & <by_project$`27`$w12> & <by_project$`30`$w12> & <by_project$`51`$w12> \\\\ \
  \\bottomrule
  \\label{tab:t1}
  \\end{tabular}
  \\end{table}",
  .open =  "<", .close = ">"
) |>
  clipr::write_clip()
