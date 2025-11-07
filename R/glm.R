glmTrainReg<-function(train_data, class_name, feature_names, ans = list(),
    debug = FALSE) {

  feature_names <- removeUniformativeFeatures(train_data, feature_names, min_examples = 1)

  fstr = paste0(class_name, "~", paste(feature_names, sep="+", collapse="+"))
  #print(fstr)
  fobj = formula(fstr)

  orig_model = suppressWarnings(glm(fobj, family = gaussian, data = train_data));

  model = orig_model;

  for (fn in feature_names) {
    if (is.factor(train_data[,fn])) {
      model$xlevels[[fn]] = levels(train_data[,fn])
    }
  }

  ans$model = model
  ans$orig_model = orig_model;
  ans$summary = summary(orig_model);
  ans$prior = mean(train_data[,class_name], na.rm=TRUE);
  ans$train_predict = suppressWarnings(predict(model, train_data, type="response"));

  return(ans);


}

glmTestReg<-function(model, test_data, feature_names, ans = list(), debug=FALSE) {
  test_p = suppressWarnings(predict(model, test_data, type="response"));
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
#'
#' @return
#' @export
#'
#' @examples
glmTrainTestReg<-function(train_data, test_data, class_name, feature_names, debug = FALSE) {
  ans = glmTrainReg(train_data, class_name, feature_names)
  ans = glmTestReg(ans$model, test_data, feature_names, ans = ans)
  return(ans);
}

