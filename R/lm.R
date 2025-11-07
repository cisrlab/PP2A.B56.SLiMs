

lmTrainReg<-function(
  train_data,
  class_name,
  feature_names,
  ans = list(),
  debug = FALSE,
  transform_fxn = "log",
  weights = rep(1, nrow(train_data)),
  var.interact = FALSE,
  include.intercept = TRUE

) {
  if (debug) {
    print(head(train_data));
    print(weights);
  }
  weights = as.vector(weights/sum(weights));
  if (debug) {
    print(weights);
  }
  #TODO, check for uninformative features and remove them.

  feature_names <- removeUniformativeFeatures(train_data, feature_names)


  fobj = lmFObj(
    class_name = class_name,
    feature_names = feature_names,
    var.interact = var.interact,
     include.intercept = include.intercept
  )
  if (debug) {
    print(fobj);
  }

  orig_model = suppressWarnings(lm(fobj, data = train_data, weights = NULL))
  model = orig_model
  for (fn in feature_names) {
    if (is.factor(train_data[, fn])) {
      model$xlevels[[fn]] = levels(train_data[, fn])
    }
  }
  ans$feature_names <- feature_names
  ans$model = model
  ans$orig_model = orig_model
  ans$summary = summary(orig_model)
  ans$train_predict = suppressWarnings(predict(orig_model, train_data, type = "response"))
  return(ans);
}

lmTestReg<-function (model, test_data, feature_names, ans = list(), debug = FALSE)
{
    test_p = suppressWarnings(predict(model, test_data, type = "response"))
    ans$test.predict = test_p
    return(ans)
}

#' Title
#'
#' @param train_data
#' @param test_data
#' @param class_name
#' @param feature_names
#' @param debug
#' @param var.interact
#' @param transform_fxn
#'
#' @return
#' @export
#'
#' @examples
lmTrainTestReg<-function(train_data, test_data, class_name, feature_names, debug = FALSE, var.interact=FALSE, transform_fxn = "none") {
    ans = lmTrainReg(train_data = train_data, class_name = class_name, feature_names = feature_names, debug=debug, var.interact = var.interact);
    ans = lmTestReg(ans$model, test_data, feature_names, ans = ans, debug=debug)
    return(ans)
}


#' Title
#'
#' @param class_name
#' @param feature_names
#' @param var.interact
#' @param include.intercept
#'
#' @return
#' @export
#'
#' @examples
lmFObj<-function(class_name, feature_names, var.interact = FALSE, include.intercept=TRUE) {
    fstr_base = paste0(class_name, "~")
    if (!include.intercept) {
      fstr_base = paste0(fstr_base, "0+")
    }


    fstr = paste0(fstr_base, paste(feature_names, sep = "+",
        collapse = "+"))
    if (var.interact) {
      int.strs = c(); #Pair-wise interactions
      for (f1 in 1:(length(feature_names)-1)) {
        for (f2 in (f1+1):length(feature_names)) {
          int.strs = c(int.strs, paste0(feature_names[f1],"*",feature_names[f2]))
        }
      }
      fstr = paste0(fstr, "+", paste(int.strs, collapse = "+", sep = "+"))
    }

  fobj = formula(fstr)
  return(fobj);
}

