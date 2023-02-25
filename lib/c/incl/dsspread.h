#ifndef DSSPREAD_HEADER
#define DSSPREAD_HEADER

/* ==== HEADER dsspread.h ==== */

/* DSSP reader routine. */

/* ANSI C, IRIX 5.3, 29. Apr. 1996. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ---- TYPEDEFS ---- */

/* Hbondrec_ : stores the offset and energy for main-chain H-bonds.
 * Residue i is bound to (i+Offs) with an energy of En kcal/mol.
 */
typedef struct
{
    int Offs;	/* 0 if unbound */
    double En;	/* 0.0 if unbound */
} Hbondrec_;

/* Dssprec_ : this struct stores the DSSP entry for a given residue. */
typedef struct 
{
    int Resno;   /* residue no. */
    char Chain;  /* chain identifier , <SPACE> if 1 chain only */
    char Res;    /* AA code, 1-letter */
    char Disulf;    /* S-S bond code or ' ' (Res set to 'C' for these) */
    char Secstruct; /* secondary structure code, see ref. */
    char Turns3, Turns4, Turns5;    /* helix turn labels */
    char Bend;	    /* geometrical bend */
    char Chir;	    /* chirality: '+' or '-' */
    char Bridge1, Bridge2;  /* beta labels: CAPS if antiparallel, lowercase if parallel */
    int Beta1, Beta2;	    /* partner resnums in beta-bridges */
    char Sheet;		    /* beta-sheet label (always caps) */
    int Access;  /* solvent accessibility, see ref. */
    Hbondrec_ Nho[2], Ohn[2];	/* H-bonds NH(i)->O and O(i)->NH */
    double Tco;	    /* CO(i):CO(i-1) torsion cosine */
    double Kappa;   /* CA(i-2):CA(i):CA(i+2) angle */
    double Alpha;   /* C-alpha dihedral angle (360 if n/a) */
    double Phi, Psi;	/* Ramachandran angles */
    double Ca[3]; /* C-alpha coordinates */
} Dssprec_;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* dssp_read: reads the text file Dsspfnm produced
 * by DSSP, Version Oct. 1988 (ref: W. Kabsch, C. Sander, Biopolymers
 * 22:2577-2637 (1983)). Returns an array of size Size, containing
 * chain ID, residue no, 1-letter AA code, Kabsch/Sander secondary
 * structure code and solvent accessibility for each residue
 * or NULL on error. Sets the number of chains to Chainno.
 */
Dssprec_ *dssp_read(const char *Dsspfnm, unsigned int *Size, unsigned int *Chainno);

/* dssp_cadist(): returns the CA:CA distance between two DSSP records
 * pointed to by Dp1 and Dp2. Checks if they're NULL.
 */
double dssp_cadist(const Dssprec_ *Dp1, const Dssprec_ *Dp2);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER dsspread.h ==== */
#endif	/* DSSPREAD_HEADER */
