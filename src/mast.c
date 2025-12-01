/************************************************************************
*                                                                       *
*       MAST                                                            *
*       Author: Timothy L. Bailey                                       *
*                                                                       *
*       Copyright                                                       *
*       (1994 - 2001) The Regents of the University of California.      *
*       All Rights Reserved.                                            *
*                                                                       *
*       Permission to use, copy, modify, and distribute any part of     *
*       this software for educational, research and non-profit purposes,*
*       without fee, and without a written agreement is hereby granted, *
*       provided that the above copyright notice, this paragraph and    *
*       the following three paragraphs appear in all copies.            *
*                                                                       *
*       Those desiring to incorporate this software into commercial     *
*       products or use for commercial purposes should contact the      *
*       Technology Transfer Office, University of California, San Diego,*
*       9500 Gilman Drive, La Jolla, California, 92093-0910,            *
*       Ph: (619) 534 5815.                                             *
*                                                                       *
*       IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO     *
*       ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR         *
*       CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF   *
*       THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA  *
*       HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
*                                                                       *
*       THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*       UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE          *
*       MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*       THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND       *
*       EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*       INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF        *
*       MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT    *
*       THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,           *
*       TRADEMARK OR OTHER RIGHTS.                                      *
************************************************************************/

/**********************************************************************/
/*
  mast <memefile> <database> [optional arguments...]
        See <install-path>/bin/mast for documentation.
*/
/**********************************************************************/

//#define DEFINE_GLOBALS
#include "mast.h"

typedef struct strand_t STRAND_T;

typedef struct sseq_t SSEQ_T;

// sortable sequence score record
struct strand_t {
    SSEQ_T *sequence; // details of the sequence
    double Pvalue; // p-value of product of p-values
    int strand; // -1 neg. strand, 0 both/protein, +1 pos. strand
    double *best_scores; // the best score for each motif
    int *best_location; // the location of the best scores for each motif
};


struct sseq_t {
    int db_index; // identifies the database
    int index; // loading number of the sequence
    long  fp; // file pointer to beginning of sequence rec.
    long length; // length of sequence
    ARRAY_T *comp; // actual sequence composition or NULL if
    int pv_alloc_len; // what is the length of the allocated pv matrix, or 0 if a reference which we don't need to free
    double **pv; // pvalue distribution for each motif (may be allocated or a reference)
    STRAND_T *pos_strand; // forward strand
    STRAND_T *neg_strand; // reverse strand
};
// }}}


/**********************************************************************/
/*
 * strand_destroy
 *
 * Destroys a strand score structure
 */
/**********************************************************************/
/*static*/ void strand_destroy(STRAND_T *strand) {
    if (strand->best_scores) myfree(strand->best_scores);
    if (strand->best_location) myfree(strand->best_location);
    memset(strand, 0, sizeof(STRAND_T));
    myfree(strand);
}


/**********************************************************************/
/*
 * sseq_destroy
 *
 * Destroys a scored sequence
 * Note, also destroys connected
 * strand score structures.
 */
/**********************************************************************/
static void sseq_destroy(void *p) {
  SSEQ_T *sseq;
  sseq = (SSEQ_T*)p;
  // deallocate connected strand objects
  if (sseq->neg_strand) strand_destroy(sseq->neg_strand);
  if (sseq->pos_strand) strand_destroy(sseq->pos_strand);
  // free up optional allocations
  if (sseq->comp) free_array(sseq->comp);
  //free expected allocations
  if (sseq->pv_alloc_len) {
    free_2array(sseq->pv, sseq->pv_alloc_len);
  }
  //zero memory
  memset(sseq, 0, sizeof(SSEQ_T));
  //free struct
  myfree(sseq);
}
