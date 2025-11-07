#' Title
#'
#' @param actual
#' @param predict
#' @param sd_actual
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
lcompplot<-function(actual,predict, sd_actual=NA, ...) {


    minx = min(min(actual[actual>0]),min(predict[predict>0]));
    maxx = max(max(actual),max(predict));

    minx = max(1e-10,minx/1.1)
    maxx = maxx * 1.1;
    #cat("minx:",minx, " ","maxx:",maxx,"\n");
    plot(actual, predict, xlim=c(minx,maxx), ylim=c(minx,maxx), log="xy", ...);


    lines(c(minx, maxx), c(minx, maxx), col="red")
    lines(c(minx, maxx), c(minx/2, maxx/2), col="purple")
    lines(c(minx, maxx), c(minx*2, maxx*2), col="purple")
    if (length(sd_actual) == length(actual)) {
        a_high = actual + 1.960 * sd_actual #95% conf assuming normal
        a_low = actual - 1.960 * sd_actual #95% conf assuming normal
        a_high[a_high <=0] = 1e-10;
        a_low[a_low <= 0] = 1e-10;
        arrows(x0 = a_low, y0 = predict, x1 = a_high, y1 = predict, col = "orange", lwd=1, angle=90, length=0.05, code=3)
        #print(cbind(actual, sd_actual, a_high, a_low))
    }

    mse <- mse(actual,predict)
    fc <- fold.err(actual,predict)

    legend("topleft", legend=c(paste("MSE", format(mse,digits=4),"FC",format(fc,digits=4))))

    ret <- mse
    attr(ret, "mse") <- mse
    attr(ret, "fc") <- fc
    invisible(ret)

}
