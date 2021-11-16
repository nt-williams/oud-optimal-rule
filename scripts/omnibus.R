library(dplyr)

source("./R/omnibus.R")
source("./R/utils.R")

# importing data
.data <- readRDS(here::here("data", "drv", "imputed-coded.rds"))[[1]] |>
  filter(medicine == "bup")

# fits <- readRDS(here::here("data", "drv", "onestep-tsm-imputed.rds")
# Y <- fits[[1]]$bup$outcome_reg[, 1]

Y <- .data$week_12_relapse
A <- .data$project
W <- .data[, w]

projects_27_30 <- .data$project %in% c("27", "30")
projects_27_51 <- .data$project %in% c("27", "51")
projects_30_51 <- .data$project %in% c("30", "51")

est.psi.prob.binom(
  W[projects_27_30, ],
  ifelse(A[projects_27_30] == "27", 1, 0),
  Y[projects_27_30],
  sig.meth = 'var',
  est.g = TRUE
) |> round(4)

est.psi.prob.binom(
  W[projects_27_51, ],
  ifelse(A[projects_27_51] == "27", 1, 0),
  Y[projects_27_51],
  sig.meth = 'var',
  est.g = TRUE
) |> round(4)

est.psi.prob.binom(
  W[projects_30_51, ],
  ifelse(A[projects_30_51] == "30", 1, 0),
  Y[projects_30_51],
  sig.meth = 'var',
  est.g = TRUE
) |> round(4)
