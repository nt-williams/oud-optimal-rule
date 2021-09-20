library(lmtp)
library(SuperLearner)

oud <- import_some_data

drug_factory <- function(x = c("methadone", "xrntx", "bupnx")) {
  x <- match.arg(x)
  \(data, trt) factor(rep(x, nrow(data)), levels = c("methadone", "xrntx", "bupnx"))
}

# Return the treatment that minimizes the conditional risk of relapse
dv <- function(risk_methadone, risk_xrntx, risk_bupnx) {
  ans <- c("methadone" = risk_methadone, "xrntx" = risk_xrntx, "bupnx" = risk_bupnx)
  names(ans[which.min(ans)])
}

# obtain the EIFs for the counterfactual risks of relapse under different trts
est_methadone <- lmtp_sdr(..., shift = drug_factory("methadone"))
est_xrntx <- lmtp_sdr(..., shift = drug_factory("xrntx"))
est_bupnx <- lmtp_sdr(..., shift = drug_factory("bupnx"))

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
