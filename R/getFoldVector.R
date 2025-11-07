#' Returns a vector that assigns the fold index for use in n-fold cross validation
#' By default do leave one out
#' @param data_in
#' @param class_name
#' @param nfold
#'
#' @return
#' @export
#'
#' @examples
getFoldVector<-function(data_in, class_name, nfold = -1) {


    if (nfold == nrow(data_in) || nfold == -1) {
        ans = sample(1:nrow(data_in), nrow(data_in));
    } else if (nfold < 1 || nfold > nrow(data_in)) {
        stop("Illegal nfold: 0 < ", nfold, "<=", nrow(data_in));
    } else {
        ans = rep(NA, nrow(data_in));

        class_values = unique(data_in[,class_name]);

        #print(class_values);

        for (class_value in class_values) {

            class_indices = which(data_in[,class_name] == class_value);
            class_folds = c();
            while ((length(class_folds)+nfold) < length(class_indices)) {
                class_folds = c(class_folds, sample(1:nfold, nfold));
            }

            rem = length(class_indices) - length(class_folds);
            if (rem) {
                class_folds = c(class_folds, sample(1:nfold, rem));
            }
            #print(class_value);
            #print(class_indices);
            #print(class_folds);
            ans[class_indices] = class_folds;

        }
    }
    if (length(ans) != nrow(data_in)) {
        print(ans);
        stop("!!!!!!");
    }

    #print(ans);
    return(ans);

}
