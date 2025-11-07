#' Given the data, class name, feature names, and the train test function of the classifier
#' Get the accuracys fo r the
#'
#' @param data_in
#' @param class_name
#' @param feature_names
#' @param train_test_fxn
#' @param nfold
#' @param fold.vec
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
getCVAccuracy<-function(
        data_in,
        class_name,
        feature_names,
        train_test_fxn = rForestTrainTest,
        nfold = -1,
        fold.vec = getFoldVector(data_in, class_name, nfold),
        ...
) {

    if (is.numeric(data_in[,class_name])) {
        #cat("CV Regression\n");
        return(getCVAccuracyFoldsReg(data_in, class_name, feature_names, train_test_fxn, fold.vec, ...));
    } else {
        return(getCVAccuracyFolds(data_in, class_name, feature_names, train_test_fxn, fold.vec, ...));
    }

}


getCVAccuracyFoldsReg<-function(data_in, class_name, feature_names, train_test_fxn, folds, track_folds=FALSE, ...) {

    ufolds = unique(folds);

    actuals = data_in[,class_name]
    predicts = rep(NA, nrow(data_in))
    ans = list();
    ans$fold.vec = folds;

    for (fold in ufolds) {
        train_data = data_in[folds != fold,];
        test_data = data_in[folds == fold,];
        fold.results = train_test_fxn(train_data, test_data, class_name, feature_names, ...);
        ans[[paste0("fold.",fold)]] = fold.results;

        predicts[folds == fold] = fold.results$test.predict;
        if (track_folds) { cat(".")};
    }
    if (track_folds) {cat("Done\n");}
    temp = predicts - actuals;
    t2 = temp * temp;

    mse = sum(t2) / length(t2);
    mae = sum(abs(temp)) / length(temp);

    nmae = mae / mean(actuals)

    perr = temp / actuals
    perr[actuals == 0] = temp[actuals == 0] / 1e-05
    perr = abs(perr) * 100;

    sstot = actuals - mean(actuals);
    sstot = sum(sstot * sstot);
    ssres = sum(t2);
    rsq = 1 - (ssres/sstot);

    ans$actuals = actuals;
    ans$predicts = predicts;
    ans$pcor = cor(predicts, actuals, method="pearson", use="complete.obs");
    ans$scor = cor(predicts, actuals, method="spearman", use="complete.obs");
    ans$mse = mse;
    ans$mae = mae;
    ans$nmae = nmae;
    ans$perr = mean(perr);
    ans$sstot = sstot;
    ans$ssres = ssres;
    ans$rsq = rsq;
    return(ans);

}


