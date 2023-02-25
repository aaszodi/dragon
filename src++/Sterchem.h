#ifndef STERCHEM_HEADER
#define STERCHEM_HEADER

// ==== PROJECT DRAGON: HEADER Sterchem.h ====

/* General stereochemical adjustment routines. */

// SGI C++ 4.0, IRIX 5.3, 22. May 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- MODULE HEADERS ----

#include "Points.h"
#include "Pieces.h"
#include "Trimat.h"

// ---- PROTOTYPES ----

/* apply_secstruct(): RMS fits all secondary structure elements in
 * Pieces onto Model, if it is 3-dimensional. 
 * Returns maximal RMS value if the average is closer to the maximum
 * than to the minimum; returns the average otherwise.
 */
double apply_secstruct(const Pieces_& Pieces, Points_& Model);

/* hand_check(): Works for 3D molecules with secondary structure only.
 * A mirror-image of the Model is generated inside and the secondary
 * structures in Pieces are fitted onto both the original and the mirror
 * image. The model with more better fits is kept, in case of a tie
 * the smaller maximal RMS value wins.
 * Return value: 1 if the original model was kept, -1 if the mirror image.
 */
int hand_check(const Pieces_& Pieces, Points_& Model);

// ==== END OF HEADER Sterchem.h ====
#endif
