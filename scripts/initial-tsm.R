library(lmtp)
library(sl3)
library(future)

# setting up parallel processing
plan(multisession)

# load functions
box::use(./R/blip, ./R/utils[...])

progressr::handlers(global = TRUE)

# importing imputed data
oud <- readRDS(here::here("data", "drv", "imputed-coded.rds"))

# creating one-step estimator function
onestep <- purrr::partial(
  lmtp_sdr,
  trt = a,
  outcome = y,
  baseline = w,
  folds = 10,
  .SL_folds = 10,
  learners_outcome = Q_learner,
  learners_trt = g_learner
)

set.seed(42523)

tsms <- list()
for (i in seq_along(oud)) {
  .data <- oud[[i]]

  # obtain the EIF for the counterfactual TSM under different medicines using one-step estimator
  met <- onestep(data = .data, shift = drug_assignment("met"))
  bup <- onestep(data = .data, shift = drug_assignment("bup"))
  nal <- onestep(data = .data, shift = drug_assignment("nal"))

  tsms[[i]] <- list(
    met = met,
    bup = bup,
    nal = nal
  )
}

# saving TSM fits and blips
saveRDS(tsms, here::here("data", "drv", "onestep-tsm-imputed.rds"))
