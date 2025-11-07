#' glmNetTrainTest
#'
#' Trains and tests a glmnet classifier
#' @param train_data training data frame with class, features as columns
#' @param test_data testing data frame with class, features as columns
#' @param class_name column name of class
#' @param feature_names vector of column names to use in the model
#'
#' @export
glmNetTrainTest <- function(train_data, test_data, class_name, feature_names, debug=FALSE) {
  library("e1071");
  ans = glmNetTrain(
    train_data = train_data,
    class_name = class_name,
    feature_names = feature_names,
    debug = debug
  );
  ans = glmNetTest(ans$model, test_data, feature_names, pans = ans, debug=debug);
  return(ans);
}

glmNetTrain<-function(
  train_data,
  class_name,
  feature_names,
  pans = list(),
  debug=FALSE) {

  ans = pans;

  fn = c();

  #glmnet doesn't like factors for some reason...
  for (feature in feature_names) {
    if (is.factor(train_data[,feature])) {
      if(nlevels(train_data[,feature]) == 2) {
        cat("Converting ",feature, " to logical\n");
        new_feature = paste0("b",feature)
        train_data[,new_feature] = train_data[,feature] == levels(train_data[,feature])[1]
        fn = c(fn, new_feature)
      } else {
        stop("Too many levels!");
      }
    } else {
      fn = c(fn, feature);
    }
  }


  class_values = train_data[,c(class_name)];

  fstr = paste0(class_name, " ~ .");
  f = formula(fstr);
  #print(f);
  library(glmnet);
  ans$model = list();
  ans$model$fit = cv.glmnet(y = train_data[,class_name], x = as.matrix(train_data[,fn]), family="binomial", type.measure="auc", nfold=3);
  plot(ans$model$fit);

  ans$train_p = predict(ans$model$fit, as.matrix(train_data[,fn]), s = "lambda.min", type="response");

  ans$roc.p = getROC_PR(data.frame(prob = ans$train_p, class=train_data[,class_name], stringsAsFactors=FALSE));
  ans$roc.n = getROC_PR(data.frame(prob = 1-ans$train_p, class=train_data[,class_name], stringsAsFactors=FALSE));

  if (getAUC(ans$roc.n)$auc_roc > getAUC(ans$roc.p)$auc_roc) {
    ans$model$reverse = TRUE;
    ans$train.prob = data.frame(true_ = 1-ans$train_p, false_ = ans$train_p, stringsAsFactors=FALSE);
  } else {
    ans$model$reverse = FALSE;
    ans$train.prob = data.frame(true_ = ans$train_p, false_ = 1-ans$train_p, stringsAsFactors=FALSE);
  }
  colnames(ans$train.prob) = c("TRUE", "FALSE");

  return(ans);

}

glmNetTest<-function(model, test_data, feature_names, pans = list(), debug = FALSE) {
  library("e1071")
  ans = pans;

  fn = c();

  for (feature in feature_names) {
    if (is.factor(test_data[,feature])) {
      if(nlevels(test_data[,feature]) == 2) {
        new_feature = paste0("b",feature)
        test_data[,new_feature] = test_data[,feature] == levels(test_data[,feature])[1]
        fn = c(fn, new_feature)
      } else {
        stop("Too many levels!");
      }
    } else {
      fn = c(fn, feature);
    }
  }




  test_p = predict(model$fit, as.matrix(test_data[,fn]), s = "lambda.min", type="response");
  if (model$reverse) {
    ans$test.prob = data.frame(true_ = 1-test_p, false_ = test_p, stringsAsFactors=FALSE);
  } else {
    ans$test.prob = data.frame(true_ = test_p, false_ = 1-test_p, stringsAsFactors=FALSE);
  }
  colnames(ans$test.prob) = c("TRUE", "FALSE");
  #print(ans$test.prob);

  return(ans);
}

glmNetGetNonZeroCoef <- function(model) {

  ans.sparse = coef(model$fit$glmnet.fit, model$fit$lambda.min)
  ans.mat = as.matrix(ans.sparse);
  nonzero = ans.mat[,1] != 0 & rownames(ans.mat) != "(Intercept)"
  nonzero.ans = ans.mat[nonzero,];
  return(nonzero.ans);

}


glmNetTrainReg<-function(train_data, class_name, feature_names, alpha=1,  ans = list(), class.wgt = "none") {
    library(glmnet)
    if (length(feature_names) == 1) {
        train_data$temp_ = rnorm(nrow(train_data), sd = 0.1);
        feature_names = c(feature_names, "temp_");
    }
    x <- glmnet::makeX(train_data[,feature_names])
    y <- train_data[,class_name]

    wgt = NULL;

    if (class.wgt == "none") {
        wgt = rep(1, nrow(train_data))
    } else if (class.wgt == "inv_linear") {
        r = nrow(train_data) - rank(train_data[,class_name]) + 1
        wgt = r;
    } else if (class.wgt == "inv_quad") {
        r = nrow(train_data) - rank(train_data[,class_name]) + 1
        wgt = r;
        wgt = wgt * wgt;
    }
    library(glmnet)
    fit <- glmnet::cv.glmnet(
        x = x,
        y = y,
        family = "gaussian",
        nfold = nrow(train_data),
        type.measure = "mse",
        parellel = TRUE,
        weights = wgt,
        grouped = FALSE,
        alpha = alpha
    )
    ans$model <- list()
    ans$model$fit <- fit
    ans$model$train_data <- train_data[,feature_names];
    ans$model$fn <- feature_names
    ans$train_predict <- predict(
        ans$model$fit,
        x,
        s = "lambda.min",
        type = "response"
    )

  return(ans);
}


#' glmnet for regression
glmNetTestReg<-function(model, test_data, feature_names, ans = list(), debug = FALSE) {
    library(glmnet)
    if (length(feature_names) == 1) {
        test_data$temp_ = rnorm(nrow(test_data), sd = 0.1)
        feature_names = c(feature_names, "temp_");
    }
    x_test <- glmnet::makeX(model$train_data, test_data[,model$fn])[["xtest"]]
    test_p <- predict(model$fit, x_test, s="lambda.min", type="response")
    ans$test.predict <- as.numeric(test_p)
    return(ans)
}


#' Title
#'
#' @param train_data
#' @param test_data
#' @param class_name
#' @param feature_names
#' @param debug
#' @param alpha
#' @param class.wgt
#'
#' @return
#' @export
#'
#' @examples
glmNetTrainTestReg<-function(train_data, test_data, class_name, feature_names, debug = FALSE, alpha=1, class.wgt = "none") {
  ans = glmNetTrainReg(train_data, class_name, feature_names, class.wgt = class.wgt, alpha = alpha);
  ans = glmNetTestReg(ans$model, test_data, feature_names, ans = ans);
}

#' Title
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
glmNetTrainTestRegRidge<-function(...) {
  glmNetTrainTestReg(alpha = 0, ...)
}




