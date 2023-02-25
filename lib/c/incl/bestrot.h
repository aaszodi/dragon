#ifndef BESTROT_H
#define BESTROT_H

/* ==== HEADER bestrot.h ==== */

/* An implementation of the point set alignment algorithm by
 * A. D. McLachlan. Reference:
 * McLachlan, A. D. (1979): J. Mol. Biol. 128: 49-79.
 * Replaces the buggy Kabsch rotation algorithm.
 */

/* ANSI C, 30. Apr. 1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ---- INCLUDE FILES ---- */

#include "matrix.h"

/* ---- PROTOTYPES ---- */

/* center_vectors: calculates the weighted centroid of the 
 * set of 3-dimensional vectors X (Vno x 3) and
 * subtracts it from each of them,  thus centring the
 * set on the centroid. If Ctr==NULL, then a 3-long
 * array is allocated to store the centroid coordinates;
 * if Ctr!=NULL, then it is assumed to be large enough to
 * hold the coordinates.
 * If W==NULL, then uniform weighting is assumed.
 * Return value: Ctr, or NULL if Vno==0.
 */
double *center_vectors(double **X, double *Ctr, const double *W,
	 unsigned int Vno);

/* best_rot: finds the best rotation matrix that brings a set of
 * vectors X into another set Y. X, Y have Vno vectors (in rows), 
 * and both live in 3 dimensions (Vno x 3). W is a Vno-long
 * weight vector that can emphasise vector pairs. 
 * If W==NULL, then uniform weighting is assumed.
 * Transform is a 3x3 square matrix (allocated before call) that on
 * return contains the X->Y transformation. It is assumed that
 * X and Y were centered before the call.
 * NOTE: this routine cannot handle the degenerate cases when
 * the point sets are Dim<3-dimensional. (Might be implemented
 * later.) When this happens, a warning is printed to stderr
 * and -1.0 (a meaningless RMS value) is returned.
 * Return value: a weighted least-squares error function.
 */
double best_rot(double **X, double **Y, const double *W, 
		unsigned int Vno, Sqmat_ Transform);

/* ==== END OF HEADER bestrot.h ==== */

#endif	/* BESTROT_H */
