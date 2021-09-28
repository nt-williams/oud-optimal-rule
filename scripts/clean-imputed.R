library(mice)
library(tidyverse)

load(here::here("data", "src", "clean_combined_imputed_data9-8-21.Rdata"))

imputed <- complete(ALT_patients_imputed_03, "all")

weeks_to_relapse <- function(rand, relapse) {
  time <- lubridate::interval(rand, relapse)
  days_to_relapse <- lubridate::as.period(time, unit = "day")
  days_to_relapse / lubridate::as.period(lubridate::dweeks(1))
}

imputed <- map(imputed, function(data) {
  mutate(
    data,
    across(starts_with("has"), ~ as.numeric(.x == "yes")),
    across(ends_with("disorder"), ~ as.numeric(.x == "yes")),
    across(c("bamphetamine30_base", "bcannabis30_base", "bbenzo30_base", "ivdrug"),
           ~ as.numeric(.x == "yes")),
    across(c("switched_meds", "never_initiated"), as.numeric),
    weeks_to_relapse = weeks_to_relapse(rand_dt, relapse_date),
    week_12_relapse = as.numeric(weeks_to_relapse <= 12),
    sex = if_else(sex == "female", 2, 1)
  )
})

saveRDS(imputed, here::here("data", "drv", "imputed-coded.rds"))
