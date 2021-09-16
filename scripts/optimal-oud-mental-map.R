library(lmtp)
library(SuperLearner)

oud <- import_some_data

# All patients recieve methadone
policy_methadone <- function(data, trt) {
  rep("methadone", length(data[[trt]]))
}

# All patients recieve XR-NTX
policy_xrntx <- function(data, trt) {
  rep("xrntx", length(data[[trt]]))
}

# All patients receive BUP-NX
policy_bupnx <- function(data, trt) {
  rep("bupnx", length(data[[trt]]))
}

# Return the treatment that minimizes the conditional risk of relapse
dv <- function(risk_methadone, risk_xrntx, risk_bupnx) {
  ans <- c("methadone" = risk_methadone, "xrntx" = risk_xrntx, "bupnx" = risk_bupnx)
  names(ans[which.min(ans)])
}

# obtain the EIFs for the counterfactual risks of relapse under different trts
est_methadone <- lmtp_tmle(..., shift = policy_methadone)
est_xrntx <- lmtp_tmle(..., shift = policy_xrntx)
est_bupnx <- lmtp_tmle(..., shift = policy_bupnx)

# regress EIFs on covariates to obtain counterfactual estimates conditional on covariates?
conditional_risk_methadone <- est_methadone$eif ~ w
conditional_risk_xrntx <- est_xrntx$eif ~ w
conditional_risk_bupnx <- est_bupnx$eif ~ w

# Apply the optimal rule to the data
oud_dv <- oud

oud_dv$trt <- dv(
  conditional_risk_methadone,
  conditional_risk_xrntx,
  conditional_risk_bupnx
)

# Estimate the risk of relapse under the optimal rule
est_dv <- lmtp_tmle(..., shifted = oud_dv)
