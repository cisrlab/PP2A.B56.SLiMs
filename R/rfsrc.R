rfsrcTrainReg<-function(train_data, class_name, feature_names, ans = list(), debug = FALSE, class.wgt = "none", seed = 1003) {
  library(randomForestSRC)
  formula_str = paste0(class_name,"~",paste(feature_names, sep="+",collapse="+"))
  #print(formula_str);
  formula_obj = as.formula(formula_str);
  fit = NULL;
  wgt = NULL;
  if (class.wgt == "none") {
    wgt = NULL
  } else if (class.wgt == "inv_linear") {
    r = nrow(train_data) - rank(train_data[,class_name]) + 1
    wgt = r
  } else if (class.wgt == "inv_quad") {
    r = nrow(train_data) - rank(train_data[,class_name]) + 1
    wgt = r
    wgt = wgt * wgt
  }

  fit <- rfsrc(
      formula_obj,
      data = train_data,
      case.wt = wgt,
      seed = -seed,
      do.trace=debug
  )

  ans$model = list();
  ans$model$fit = fit;
  ans$model$wgt = wgt;
  ans$model$formula_str = formula_str;
  ans$model$formula_obj = formula_obj;
  ans$train_predict = predict(fit, train_data)$predicted;
  return(ans);
}

#' Title
#'
#' @param model
#' @param test_data
#' @param feature_names
#' @param ans
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
rfsrcTestReg<-function(model, test_data, feature_names, ans = list(), debug = FALSE) {
  test_p = predict(model$fit, test_data)$predicted;
  ans$test.predict = test_p;
  return(ans);
}


#' Title
#'
#' @param train_data
#' @param test_data
#' @param class_name
#' @param feature_names
#' @param debug
#' @param class.wgt
#'
#' @return
#' @export
#'
#' @examples
rfsrcTrainTestReg<-function(train_data, test_data, class_name, feature_names, debug = FALSE, class.wgt = "none") {
  ans = rfsrcTrainReg(train_data, class_name, feature_names, class.wgt = class.wgt);
  ans = rfsrcTestReg(ans$model, test_data, feature_names, ans = ans);
  return(ans);
}

rfsrcTrainTestW1Reg<-function(train_data, test_data, class_name, feature_names, debug=FALSE) {
  return(rfsrcTrainTestReg(train_data, test_data, class_name, feature_names, debug = debug, class.wgt = "inv_linear"))
}

rfsrcTrainTestW2Reg<-function(train_data, test_data, class_name, feature_names, debug=FALSE) {
  return(rfsrcTrainTestReg(train_data, test_data, class_name, feature_names, debug = debug, class.wgt = "inv_quad"))
}


rfsrc.min.mperr<-function(data_in, wgt) {
  forumula_obj = attr(data_in, "formula_obj")
  actual_col = attr(data_in, "actual_col")
  transform_fxn = attr(data_in, "transform_alg")
  seed = attr(data_in, "seed")
  square = attr(data_in, "square")
  fit = try(rfsrc(formula= forumula_obj, data = data_in, case.wt = wgt));
  if ("try-error" %in% class(fit)) {
    return (NA);
  }
  predicted = predict(fit, data_in)$predicted;
  predicted = undoTransform(predicted, transform_fxn)
  actual = undoTransform(data_in[,actual_col], transform_fxn)


  train_mperr = mperr(actual = actual, predicted = predicted, square = square);
  return(train_mperr);

}

rfsrc.min.err<-function(data_in, wgt) {
  forumula_obj = attr(data_in, "formula_obj")
  actual_col = attr(data_in, "actual_col")
  transform_fxn = attr(data_in, "transform_alg")
  seed = attr(data_in, "seed")
  err_fxn = attr(data_in, "err_fxn")

  train_err = 1e10;

  fit = try(rfsrc(formula= forumula_obj, data = data_in, case.wt = wgt, seed = -seed));

  if (inherits(fit, "try-error")) {
    train_err = 1e10;
  } else {
    predicted = predict(fit, data_in)$predicted;
    predicted = undoTransform(predicted, transform_fxn)
    actual = undoTransform(data_in[,actual_col], transform_fxn)
    train_err = err_fxn(actual = actual, predicted = predicted);
  }
  return(train_err);
}


rfsrcTrainWEReg = function (train_data, class_name, feature_names, ans = list(),
    debug = FALSE, iter=10, transform_alg="none", err_fxn = mperr, seed = 1003) {
    library(randomForestSRC)
    #Build the formula object
    formula_obj <- lmFObj(class_name, feature_names)

    #Set up the parameters and data for optimization
    wgt = rep(1, nrow(train_data))
    wgt = wgt / sum(wgt);
    attr(train_data, "formula_obj") <- formula_obj;
    attr(train_data, "actual_col") <- class_name;
    attr(train_data, "transform_alg") <- transform_alg;
    attr(train_data, "seed") <- seed;
    attr(train_data, "err_fxn") <- err_fxn;

    #Get the initial error with equal weights
    err0 <- rfsrc.min.err(train_data, wgt)

    #Try to optimize weights to reduce error
    lb = rep(1e-5, nrow(train_data))
    ub = rep(1,nrow(train_data))
    res = optim(
        par = wgt,
        fn = rfsrc.min.err,
        data = train_data,
        control=list(maxit=iter,trace=as.integer(debug)),
        method = "L-BFGS-B", lower=lb, upper=ub)
    #print(res);
    #Determine best weights including the initial equal weights
    wgt_best <- wgt
    err_best <- err0

    if (res$convergence == 1) {
        if (debug) {message("maxit reached:", iter)}
    }

    if (res$value < err0) {
        wgt_best <- res$par
        err_best <- res$value
    }

    message("err0:", err0, " best:", err_best)

    #Fit the final model with the optimized weights
    fit_best <- rfsrc(
        formula_obj, data = train_data, case.wt = wgt_best, seed = -seed)

    #message("Best mperr:", mperr_best);
    ans$model = list();
    ans$model$fit = fit_best
    ans$model$err_best = err_best
    ans$model$wgt_best = wgt_best
    ans$model$formula_obj <- formula_obj
    ans$train_predict <- predict(fit_best, train_data)$predicted
    return(ans)
}




rfsrcTrainWMPReg = function (train_data, class_name, feature_names, ans = list(),
    debug = FALSE, iter=10, transform_alg="none", seed = 1003) {
    library(randomForestSRC)

    return(
      rfsrcTrainWEReg(
        train_data = train_data,
        class_name = class_name,
        feature_names = feature_names,
        ans = ans,
        debug = debug,
        iter = iter,
        transform_alg = transform_alg,
        seed = seed,
        err_fxn = mperr
      )
    )
}

rfsrcTrainTestWMPReg <- function (train_data, test_data, class_name, feature_names, transform_alg = "none", debug = FALSE) {
  ans = rfsrcTrainWMPReg(train_data, class_name, feature_names, transform_alg = transform_alg);
  ans = rfsrcTestReg(ans$model, test_data, feature_names, ans = ans)
  return(ans);
}

rfsrcTrainTestWEReg<-function(
    train_data,
    test_data,
    class_name,
    feature_names,
    transform_alg = "none",
    debug = FALSE,
    err_fxn = mperr) {

  ans = rfsrcTrainWEReg(
    train_data = train_data,
    class_name = class_name,
    feature_names = feature_names,
    transform_alg = transform_alg,
    debug = debug,
    err_fxn = err_fxn
  )
  ans = rfsrcTestReg(ans$model, test_data, feature_names, ans=ans)
  return(ans)
}


