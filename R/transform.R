#' Title
#'
#' @param x
#' @param alg
#'
#' @return
#' @export
#'
#' @examples
doTransform<-function(x, alg="none") {
  ans = x;
  if (alg == "none") {
    ans = x;
  } else if (alg =="log") {
    ans = log2(x)
  } else if (alg == "nlog") {
    ans = -log2(x)
  } else if (alg == "log1p") {
    ans = log1p(x)
  } else if (alg == "nlog1p") {
    ans = -log1p(x)
  } else if (alg == "sqrt") {
    ans = sqrt(x);
  } else if (alg == "inv") {
    ans = 1/x;
  } else {
    stop("Unknown transform:", alg);
  }
#  print(ans);
  return(ans);
}

#' Title
#'
#' @param x
#' @param alg
#'
#' @return
#' @export
#'
#' @examples
undoTransform<-function(x, alg="none") {
  ans = x;
  if (alg == "none") {
    ans = x;
  } else if (alg == "log") {
    ans = 2^x;
  } else if (alg == "nlog") {
    ans = 2^-x;
  } else if (alg == "log1p") {
    ans = expm1(x);
  } else if (alg == "nlog1p") {
    ans = expm1(-x);
  } else if (alg == "sqrt") {
    ans = x^2
  } else {
    stop("Unknown transform:", alg);
  }
  return(ans);
}

