#' Title
#'
#' @param seq
#' @param meme_xml
#' @returns
#' @export
#'
#' @examples
fimo <- function(fimo_args) {
    co <- capture.output({
    ret<- .Call("C_fimo", fimo_args)
    })
    attr(co,"ret") <- ret
    return(co)
}


#' Title
#'
#' @param fasta
#' @param meme_xml
#' @param bg_file
#' @param dir
#' @param no_qvalue
#' @param thresh
#' @param max_stored_scores
#' @param debug
#'
#' @returns
#' @export
#'
#' @examples
fimoF<-function(
    fasta,
    meme_xml,
    bg_file="",
    dir = "fimo_output",
    no_qvalue = TRUE,
    thresh=0.5,
    max_stored_scores = 100000,
    debug = FALSE
) {

  if (debug) {cat("meme$fimoFasta: start\n");}

  #fimo_cmd = meme$getFimo();

  fimo_args <- c()

  fimo_args <- c(fimo_args, "--norc")
  fimo_args <- c(fimo_args, "--thresh", thresh)
  fimo_args <- c(fimo_args, "--oc")
  if (get_os() == "windows") {
    fimo_args <- c(fimo_args, cygwinpath(dir))
  } else {
    fimo_args <- c(fimo_args, paste0("\"", dir, "\""))
  }
  if (bg_file != "") {
    fimo_args <- c(fimo_args, "--bgfile")
    fimo_args <- c(fimo_args, bg_file)
  }
  if (no_qvalue) {
    fimo_args <- c(fimo_args, "--no-qvalue")
  }
  fimo_args <- c(fimo_args, "--max-stored-scores", max_stored_scores)
  fimo_args <- c(fimo_args, "--text")

  if (get_os() == "windows") {
    fimo_args = c(fimo_args, cygwinpath(meme_xml), cygwinpath(fasta))
  } else {
    fimo_args = c(fimo_args, meme_xml, fasta)
  }

  fimo_args <- c(fimo_args, "--verbosity", 1) # QUIET_VERBOSE = 1

  #print(fimo_args)
  ret <- fimo(fimo_args)
  ans <- readFimoLines(ret)
  attr(ans, "ret") <- attr(ret, "ret")
  return(ans)
}

