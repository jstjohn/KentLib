/* 
 * File:   pssm.h
 * Author: TCO
 *
 * Created on 3. Februar 2009, 13:14
 */

#ifndef _PSSM_H
#define	_PSSM_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "params.h"
#include "string.h"

/* Define PSSM as an array of position specific substitution
   matrices to be used at different points in an alignment.
   The first PSSM_DEPTH-1 matrices are for the beginning of the
   alignment. The matrix at PSSM_DEPTH is for the middle. The
   matrices at PSSM_DEPTH+1..2*PSSM_DEPTH are for the end of
   the alignment.
   The second value is for the reference base. The third value
   is for the ancient base.
   A=0, C=1, G=2, T=3, anthing else=4
*/
typedef struct pssm {
  int sm[2*PSSM_DEPTH+1][5][5];
  int depth;
} PSSM;
typedef struct pssm* PSSMP;

    /* revcom_submat
       Takes a PSSMP (sm) pointer to a valid submat
       Makes a reverse complement of this submat
       Returns a pointer to the new submat
     */
    PSSMP revcom_submat(PSSMP psm);

    /* Reads in this hardcoded flat substitution matrix */
    PSSMP init_flatsubmat(void);

    /* s1b is the reference base
       s2b is the fragment (ancient) base */
    int sub_mat_score(const short int s1i,
            const short int s2i,
            int sm[][5][5],
            const int row,
            const int len);


    /* find_sm_depth
       Args: (1) int row - the current row we're on for alignment, i.e.
                           the position in the fragment sequence
             (2) len len - the length of the fragment sequence we are
                           aligning
       Returns: int - the depth in the substitution matrix for this position
                in the fragment
     */
    int find_sm_depth(int row, int len);

    PSSMP read_pssm( const char* fn );

#ifdef	__cplusplus
}
#endif

#endif	/* _PSSM_H */

