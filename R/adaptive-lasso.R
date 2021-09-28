.data <- as.data.frame(data_cat_realistic)[, w]
.ndata <- as.data.frame(model.matrix(~ . ^ 2, .data)[, -1])
.ndata <- cbind(.ndata, as.data.frame(data_cat_realistic)[, paste0("blip", 1:2)])

res <- lm(cbind(blip1, blip2) ~ ., data = .ndata)
pf <- (1 / abs(res$coefficients))

mstar <- glmnet::cv.glmnet(as.matrix(.ndata[, 1:10]),
                  as.matrix(.ndata[, paste0("blip", 1:2)]),
                  family = "mgaussian",
                  penalty.factor = pf[-1, 1])

rownames(coef(mstarfit, s="lambda.min"))[coef(mstarfit, s="lambda.min")[,1]!=0]

estimate_blip_multi_lasso <- function(data, covar, blip, nfolds) {
  folds <- origami::make_folds(data, V = nfolds)
}

fit_blip_multi_lasso <- function(fold, data, covar, blip) {
  train <- origami::training(data, fold)
  valid <- origami::validation(data, fold)

  browser()

  train_y <- train[, blip]
  train_x <- train[, covar]

  init <- lm.fit(train_x, train_y)
}
