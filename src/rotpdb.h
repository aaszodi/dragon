#ifndef ROTPDB_HEADER
#define ROTPDB_HEADER

/* ==== HEADER rotpdb.h ==== */

/* Contains routines for PDB I/O and for optimal superposition
 * of PDB-derived C-alpha chains.
 */

/* ANSI C, 7. Feb. 1997. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ---- MODULE HEADERS ---- */

#include "matrix.h"	/* matrix operations */
#include "pdbprot.h"   /* PDB data and I/O */
#include "bestrot.h"   /* rigid body rotation (McLachlan-algorithm) */

/* ---- PROTOTYPES ---- */

/* get_vectors: transfers the coordinates of the
 * first chain in the PDB file Pdbfnm to an array of 3D
 * vectors suitable for processing by the Bestrot routines.
 * If *Newentry==NULL, then it will be set to point to
 * the PDB entry record constructed from the contents of
 * the PDB file (used for reading the first entry). If *Newentry!=NULL,
 * then the first chain in it is compared to the first chain
 * just read in and NULL is returned if they don't match.
 * Ignored if Newentry==NULL.
 * If Allatoms!=0, then all atoms are read, otherwise C-alphas only.
 * The vectors are centred on the centroid.
 * The no. of vectors is returned in Size, the return value
 * is the [Size][3] array of the vectors or NULL if the file
 * could not be read or was not a protein.....etc.
 * Also returns the sequence of the chain in 1-letter codes in Seq.
 * Also creates and returns a weight vector in *Wgt if Wgt!=NULL;
 * the weights will be "reciprocal" to the B-factors in the
 * structure (if meaningful B-factors are specified) or
 * uniform otherwise.
 */
double **get_vectors(const char *Pdbfnm, Pdbentry_ **Newentry, 
	char Allatoms, double **Wgt, int *Size);

/* rotate_vectors: performs the McLachlan rotation on the vector sets
 * Target and Struct (both [Size][3]) so that Struct will be
 * rotated to Target. Weights are supplied in W[].
 * Return value: the RMS difference.
 */
double rotate_vectors(double **Target, double **Struct, 
	const double W[], int Size);

/* smooth_chains: do a moving-average smoothing of the
 * coordinate vectors stored in Struct (Size x 3) Cycno times
 * with a Wlen-wide window. The smoothed results will be kept in
 * Struct which will be truncated appropriately on exit.
 * Return value: a pointer to the realloc()-d Struct whose
 * new size is adjusted in Size.
 */
double **smooth_chains(double **Struct, int *Size, int Wlen, int Cycno);

/* smooth_wgt: smooths the weight vector Wgt. The algorithm and
 * the parameters used are the same as in smooth_chains().
 * Note that Size is not modified and Wgt is not-reallocated.
 */
void smooth_wgt(double Wgt[], int Size, int Wlen, int Cycno);

/* start_struct: creates a PDB entry structure (cf. "pdbprot.h")
 * from the Target vectors interpreted as coordinates
 * for the first chain in Pdbtarg. If it is NULL,
 * it means that smoothing on C-alphas was
 * done and no sequences are saved. (Residue IDs will be "XXX".)
 * Return value: ptr to the entry or NULL if Targetsize<=0.
 */
Pdbentry_ *start_struct(double **Target, int Targetsize, 
	const Pdbentry_ *Pdbtarg);

/* add_struct: adds a new chain to Entry which contains the 
 * vectors in Struct,  interpreted as coordinates belonging
 * to the first chain in Pdbtarg if it is not NULL.
 * The new chain is Targetsize long, with Chainid as a
 * character ID. It is assumed that there's a 1:1 correspondence
 * between the vectors in Struct and the C-alpha coordinates
 * in Target (see above) stored in Entry->Chain[0].
 * Structseq==NULL means that smoothing has been done and no
 * sequence information will be saved.
 * The B-factor field will contain the distance (in angstroms)
 * between the corresponding atoms in Target and Struct.
 * This is for the Quanta "4th-parameter" colouring option.
 */
void add_struct(Pdbentry_ *Entry, double **Struct, int Targetsize, 
	const Pdbentry_ *Pdbtarg, char Chainid);

/* target_sd(): calculate the standard deviation of distances 
 * from the target for each atom and put these into the corresponding
 * B-factor entries of the target chain. We assume that all chains
 * have already been added to Entry before the call and the structure atom
 * B-factor fields contain the distance from the target atom. The B-factor
 * values are set to 0.0 if there were less than 2 structures in
 * addition to the target.
 */
void target_sd(Pdbentry_ *Entry);

/* ==== END OF HEADER rotpdb.h ==== */

#endif	    /* ROTPDB_HEADER */
