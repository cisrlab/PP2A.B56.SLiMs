glmNetTrainReg2<-function(train_data, class_name, feature_names, ans = list(), class.wgt = "none") {

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
    library(glmnetUtils)

    tune_res <<- glmnetUtils::cva.glmnet(
        x = x,
        y = y,
        alpha = seq(0, 1, len = 21)^3,
        family = "gaussian",
        nfold = nrow(train_data),
        type.measure = "mse",
        parallel = TRUE,
        weights = wgt,
        grouped = FALSE
    )
    #print(tune_res)
    #stop()
    alpha_min_cvm <- unlist(lapply(tune_res$modlist,function(l) {min(l$cvm)}))
    opt_alpha_idx <- which(alpha_min_cvm == min(alpha_min_cvm))[1]


    ans$model <- list()

    ans$x <- x
    ans$y <- y
    ans$model$tune_res <- tune_res
    ans$model$alpha_min_cvm <- alpha_min_cvm
    ans$model$opt_alpha_idx <- opt_alpha_idx
    ans$model$opt_alpha <- tune_res$alpha[opt_alpha_idx]
    ans$model$fit <- tune_res$modlist[[opt_alpha_idx]]
    ans$model$train_data <- train_data[,feature_names]
    ans$train_predict <- predict(
        ans$model$fit,
        x,
        s = "lambda.min",
        type = "response"
    )

    return(ans)
}

glmNetTestReg2 <- function(model, test_data, feature_names, ans = list(), debug = FALSE) {
    library(glmnet)
    if (length(feature_names) == 1) {
        test_data$temp_ = rnorm(nrow(test_data), sd = 0.1);
        feature_names = c(feature_names, "temp_");
    }
    x_test <- makeX(model$train_data, test_data[,feature_names])[["xtest"]]
    test_p = predict(model$fit, x_test, s = "lambda.min",
                     type = "response");
    ans$test.predict = as.numeric(test_p);
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
glmNetTrainTestReg2 <-
    function(train_data, test_data, class_name, feature_names, debug = FALSE, class.wgt = "none") {
        ans <- glmNetTrainReg2(train_data, class_name, feature_names, class.wgt = class.wgt);
        ans <- glmNetTestReg2(ans$model, test_data, feature_names, ans = ans);
        return(ans)
    }
