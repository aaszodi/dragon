#ifndef __HMOM_HEADER__
#define __HMOM_HEADER__

// ==== PROJECT DRAGON: HEADER Hmom.h ====

/* Calculates the hydrophobic moments of clusters and
 * rotates them in Euclidean space so that the moments point
 * towards the common centroid.
 */

// SGI C++ 4.0, IRIX 5.3, 3. July 1995. Andris Aszodi 

// ---- UTILITY HEADERS ----

#include "Points.h"

// ---- MODULE HEADERS ----

#include "Polymer.h"
#include "Pieces.h"

// ---- PROTOTYPES ----

/* hmom_clurot(): rotates all clusters in Pieces so that their
 * hydrophobic moments point towards the overall centroid.
 * The phobicity info is in Polymer, the Euclidean coordinates
 * are in Xyz. No action is taken for single-cluster sets.
 */
void hmom_clurot(const Pieces_& Pieces, const Polymer_& Polymer, 
	Points_& Xyz);

// ==== END OF HEADER Hmom.h ====
#endif
