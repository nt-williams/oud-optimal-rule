library(lmtp, lib = "~/Desktop/lmtp-alternates/sl3")
library(sl3)
library(future)

# load functions
box::use(../R/blip, ../R/utils[...])

progressr::handlers(global = TRUE)

# importing imputed data
oud <- readRDS(here::here("data", "drv", "imputed-no-27bup.rds"))

oud |>
  (\(x) x[[1]])() |>
  (\(x) x$project == "27")() -> in27

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

plan(multisession)

tsms <- list()
for (i in seq_along(oud)) {
  .data <- oud[[i]][!in27, ]

  # obtain the EIF for the counterfactual TSM under different medicines using one-step estimator
  bup <- onestep(data = .data, shift = drug_assignment("bup"))
  nal <- onestep(data = .data, shift = drug_assignment("nal"))

  tsms[[i]] <- list(
    bup = bup,
    nal = nal
  )
}

plan(sequential)

# saving TSM fits
saveRDS(tsms, here::here("data", "drv", "onestep-tsm-imputed-no27.rds"))
