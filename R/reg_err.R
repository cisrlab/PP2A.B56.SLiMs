
#' Title
#'
#' @param actual actual numeric values
#' @param predicted predicted numeric values
#' @param square square the error?
#'
#' @return
#' @export
#'
#' @examples
mperr<-function(actual, predicted, square=FALSE) {
  temp = (predicted - actual)
  perr = temp/actual
  perr[actual == 0] = temp[actual == 0]/1e-05
  perr = abs(perr);
  perr = perr * 100;
  if (square) {
    perr = perr ^ 2
  }
  ans = mean(perr, na.rm=TRUE);
  attr(ans,"perr") = perr;
  return(ans);
}

#' Title
#'
#' @param actual
#' @param predicted
#'
#' @return
#' @export
#'
#' @examples
mse <-function(actual, predicted) {
  ans = actual - predicted
  ans = ans * ans;

  mans <- mean(ans, na.rm=TRUE)
  attr(mans, "mse") <- ans
  return(mans)


}

#' Title
#'
#' @param actual
#' @param predicted
#'
#' @return
#' @export
#'
#' @examples
mae<-function(actual, predicted) {
  ans = abs(actual - predicted)
  return(mean(ans))
}

#' Title
#'
#' @param actual
#' @param predicted
#'
#' @return
#' @export
#'
#' @examples
rmse<-function(actual, predicted) {
  ans = sqrt(mse(actual, predicted))
  return(ans);
}

#' Title
#'
#' @param actual
#' @param predicted
#'
#' @return
#' @export
#'
#' @examples
nmae<-function(actual, predicted) {
  ans = mae(actual, predicted)/ mean(actual);
  return(ans);
}

#' Mean count of the number of points whose ratio of
#'
#' @param actual
#' @param predicted
#' @param pthreshold
#'
#' actual to predict is > ratio, like a tube.
#' ratio of 2 means within 50 percent error.
#'
#' @return
#' @export
tube.err<-function(actual, predicted, pthreshold=50) {

  mperr_res = mperr(actual,predicted)
  perr_res = attr(mperr_res, "perr")
  ans = (1-mean(perr_res <= pthreshold))*100;
  attr(ans, "mperr") = mperr_res;
  return(ans);
}

tube.err2<-function(actual, predicted, pthreshold = 50) {
  mperr_res = mperr(actual,predicted)
  perr_res = attr(mperr_res, "perr")

  ans = perr_res / pthreshold;
  ans[perr_res <= pthreshold] = 0;
  ans.mean = mean(ans) * 100;
  attr(ans.mean, "mperr") = mperr_res;
  attr(ans.mean, "tube.err") = ans;
  return(ans.mean);

}


#' Title
#'
#' @param actual
#' @param predicted
#'
#' @return
#' @export
#'
#' @examples
med.perr<-function(actual, predicted) {
  temp = (predicted - actual)
  perr = temp/actual
  perr[actual == 0] = temp[actual == 0]/1e-05
  perr = abs(perr);
  perr = perr * 100;
  ans = median(perr, na.rm=TRUE);
  return(ans);

}



#' Title
#'
#' @param actual
#' @param predicted
#' @param actual_sd
#' @param square
#'
#' @return
#' @export
#'
#' @examples
avg.zerr<-function(actual, predicted, actual_sd, square=FALSE) {
  zerr = abs(predicted-actual) / actual_sd;
  if (square) {
    zerr = zerr * zerr
  }
  ans = mean(zerr, na.rm=TRUE);
  attr(ans, "zerr") = zerr;
  return(ans);
}

#' fold-change error
#'
#' @param actual
#' @param predicted
#'
#' @return
#' @export
#'
#' @examples
fold.err<-function(actual, predicted) {
    fc = pmax(actual,predicted) / pmin(actual, predicted)
    ans = mean(fc, na.rm=TRUE);
    attr(ans, "fc") = fc;
    return(ans);
}

rsq <- function(actual, predicted) {
    ma <- mean(actual, na.rm=TRUE)
    adiff <- actual - ma
    sstot <- sum(adiff * adiff)

    ediff <- predicted - actual
    sserr <- sum(ediff * ediff)

    ans <- min(max(0, 1 - sserr / sstot), 1.0)
    attr(ans, "sserr") <- sserr
    attr(ans, "sstot") <- sstot
    attr(ans, "ssreg") <- sstot - sserr
    return(ans)
}


