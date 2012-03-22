/*
 * params.h
 *
 *  Created on: Feb 28, 2012
 *      Author: jstjohn
 */

#ifndef PARAMS_H_
#define PARAMS_H_
#include <stdio.h>
#include <stdlib.h>

#define PSSM_DEPTH (15)
#define GOP (1000) // Gap open penalty
#define GEP (200) // Gap extension penalty
#define N_SCORE (-100)
#define NR_SCORE (-10) // score for N in reference
#define FLAT_MATCH (200) // score in the flat matrix for a gap
#define FLAT_MISMATCH (-600) // score in the flat matrix for a mismatch
#define MAX_LINE_LEN (1000000)
#define MAX_NAME_LEN (50)

#endif /* PARAMS_H_ */
