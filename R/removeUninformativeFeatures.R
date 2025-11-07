#' Title
#'
#' @param fs
#' @param feature_names
#'
#' @return
#' @export
#'
#' @examples
removeUniformativeFeatures<-function(
        fs,
        feature_names,
        min_examples = 1,
        debug=TRUE,
        na.rm = FALSE
) {
    to_remove = c();
    to_keep = c();
    for (feature in feature_names) {
        vals <- fs[,feature]
        if (na.rm) {vals <- vals[!is.na(vals)]}
        #Remove all features that do not have at least a FALSE and TRUE
        if (is.logical(fs[,feature])) {
            tbl <- table(vals)
            tbl <- tbl[tbl>= min_examples]
            if (length(tbl)<2) {
                to_remove <- c(to_remove, feature)
            } else {
                to_keep <- c(to_keep, feature)
            }
        } else if (is.factor(fs[,feature]) || is.character(fs[,feature])) {
            #Remove all features that do not have at least two factors with
            #Some data.
            tbl <- table(vals)
            tbl <- tbl[tbl >= min_examples]
            if (length(tbl)< nlevels(fs[,feature])) { #
                to_remove <- c(to_remove, feature)
                if (debug) {
                    print(feature)
                    print(table(fs[,feature]))
                }
            } else {
                to_keep <- c(to_keep, feature)
            }

        } else if (is.numeric(fs[,feature])) {
            #Remove all features that have only 1 numeric value
            tbl <- table(vals)
            tbl <- tbl[tbl >= min_examples]
            if (length(unique(fs[,feature]))==1) {
                to_remove <- c(to_remove, feature)
            } else {
                to_keep <- c(to_keep, feature)
            }
        } else {
            stop("Feature is not a numeric, factor or logical:", feature)
        }
    }

    stopifnot(length(to_remove)+length(to_keep) == length(feature_names))
    if (debug && length(to_remove) > 0) {
        message("Out of ", length(feature_names))
        message("Removed ", length(to_remove), " uninformative features");
        message("Kept", length(to_keep), " features")
    }

    fn <- to_keep
    attr(fn, "to_remove") = to_remove
    attr(fn, "to_keep") <- to_keep
    return(fn)
}
