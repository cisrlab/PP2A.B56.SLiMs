#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "fimo.h"

#include <stdio.h>

/*   seqs,
 meme_xml,
 bg_file="",
 dir=paste0(getwd(),"/fimo_output"),
 no_qvalue = FALSE,
 thresh=0.5,
 max_stored_scores = 100000,
 debug = FALSE
*/

SEXP C_fimo(
  SEXP fimo_args
  ) {
    // Basic type check
    if (TYPEOF(fimo_args) != STRSXP) {
        Rf_error("Expecting a character vector (STRSXP) for fimo_args");
    }
    int len = Rf_length(fimo_args);

    char **argv = malloc(len * sizeof(char *));


    for (int i = 0; i < len; i++) {
        argv[i] = CHAR(STRING_ELT(fimo_args, i));
    }


    int ret = fimo_main(len, argv);

    SEXP ret_val;
    PROTECT(ret_val = NEW_INTEGER(1));

    // Get a pointer to the underlying C integer array and set the value
    INTEGER(ret_val)[0] = ret;

    // Unprotect the SEXP before returning
    UNPROTECT(1);
    return ret_val;

}
