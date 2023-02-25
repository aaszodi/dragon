#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

// ==== PROJECT DRAGON: HEADER Output.h ====

/* Lists the simulation results to a file in PDB format. 
 * Converts the result first to the C structure Pdbentry_
 * and then uses the C module "pdbprot" for the actual output.
 */

// SGI C++ 4.0, IRIX 5.3, 10. July 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- UTILITY HEADERS ----

#include "String.h" 

// ---- MODULE HEADERS ----

#include "Polymer.h"
#include "Pieces.h"
#include "Score.h"

// ---- PROTOTYPES ----

/* make_outname(): constructs an output filename of the form
 * "Basename_X.Ext" where X is the run number Rcyc. 
 * If Basename contains a directory path, then this path
 * will be created if necessary (permissions permitting).
 * If path creation fails, then the dirpath is thrown away
 * from Basename (will be used in the current directory).
 */
void make_outname(String_& Basename, int Rcyc, const String_& Ext);

/* pdb_result(): saves the result of the simulation in a file Pdbf
 * provided the coordinates in Xyz are 3-dimensional and there are
 * no dimension mismatches. Returns 1 on success, 0 on error.
 */
int pdb_result(const char *Pdbf, const Points_& Xyz, 
	const Polymer_& Model, const Pieces_& Pieces, 
	const Scores_& Bestsco);

// ==== END OF HEADER Output.h ====

#endif	/* OUTPUT_HEADER */
