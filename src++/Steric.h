#ifndef STERIC_CLASS
#define STERIC_CLASS

// ==== PROJECT DRAGON: HEADER Steric.h ====

/* Steric adjustment routines. */

// SGI C++ 7.1, IRIX 6.2, 1. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- UTILITY HEADERS ----

#include "Trimat.h"
#include "Points.h"

// ---- MODULE HEADERS ----

#include "Polymer.h"
#include "Fakebeta.h"
#include "Restr.h"
#include "Pieces.h"
#include "Specgrad.h"
#include "Score.h"
#include "Viol.h"

// ==== CLASSES ====

/* Steric_: a class which performs most steric adjustments.
 * "Ideal" distances are calculated so that they fall within
 * distance bounds. In DIST space, these are applied to the
 * distance matrix directly, in EUCL space, through majorization
 * algorithms.
 */
class Steric_
{
    // data
    private:
    
    Trimat_ Strimat;	// the strictness matrix (weights)
    Trimat_ Idist;    // UNsquared ideal distances
    Specgrad_ Sp;   // spectral gradient for the whole lot
    int Lastflags;  // the last adjustments
    
    public:
    /* The actions of the adjustment routines are controlled by the
     * following flags:-
     * WITHIN: adjustment is done inside the clusters only;
     * BETWEEN: adjustment is done between the clusters only;
     * ALL: adjustment is done everywhere (== WITHIN | BETWEEN)
     * SPECGRAD: adjustment to be done with the Spectral Gradient refinement
     * RINT: do the bond and bump (internal) restraints only
     * REXT: do the external restraints only
     * RESTR: do all restraints (==RINT | REXT)
     * BOND: do the virtual CA:CA bonds and CA(i):CA(i+2) only
     * SCORE: generate scores
     * The flags may be OR-ed together.
     */
    enum Adjflags_ {WITHIN=1, BETWEEN=2, ALL=3, SPECGRAD=4, 
	    RINT=8, REXT=16, RESTR=24, BOND=32, SCORE=64};
    
    // methods
    public:
    
	// Constructors
    /* Inits for Resno residues (default 10). 
     * There are 2 extra points for the N/C-terminal moieties.
     */
    Steric_(unsigned int Resno=10): 
	Strimat(Resno+2), Idist(Resno+2), Lastflags(0) {}
    
	// Setup
    /* setup(): Changes the size of the matrices so that
     * they can hold data for Rno residues (Rno+2 points with the
     * N/C termini).
     * Returns old residue no. 
     */
    unsigned int setup(unsigned int Rno)
    {
	unsigned int Oldrno=Idist.rno()-2;	// old size
	Strimat.set_size(Rno+2);
	Idist.set_size(Rno+2);
	return(Oldrno);
    }
    
	// Ideal distances
    /* ideal_dist(): fills up the ideal distance matrix within the
     * calling object. Dista is the actual CA:CA distance matrix (squared), 
     * Fakebeta can be queried for sidechain centroid (SCC) distances, 
     * Restraints holds the external restraints and various bump lengths, 
     * Polymer supplies CA:SCC distances, Pieces provides the cluster
     * layout, Checkflags controls the adjustment.
     * If the SCORE flag was specified in Checkflags, then *Scores will contain
     * a score update on return, otherwise it is left unchanged.
     * If SCORE is set and Vl!=NULL then the violation list pointed to by Vl is made.
     */
    void ideal_dist(const Trimat_& Dista, const Fakebeta_& Fakebeta, 
	    const Restraints_& Restraints, const Polymer_& Polymer, 
	    const Pieces_& Pieces, int Checkflags, Scores_* Scores=NULL, Viollist_ *Vl=NULL);

	// Violation assessment
    /* reset_viol(): sets the score normalisation factors (the appropriate
     * sums of various restraint weights) in Scores.
     * Call once before a simulation run.
     */
    void reset_viol(const Restraints_& Restraints, int Size, Scores_& Scores) const;
    
	// Adjustments
    /* adjust_dist(): adjusts steric clashes in 'dist' space. Replaces the
     * entries in Dista by the corresponding entries in Idist
     * according to the cluster structure and the check choice.
     * The extent of the adjustment is controlled by Strict:
     * 0.0 means no adj, 1.0 total adj (default). 
     */
    void adjust_dist(Trimat_& Dista, const Pieces_& Pieces, 
	int Checkflags, float Strict=1.0) const;
    
    /* adjust_xyz(): tries to manoeuvre C-alphas to positions in Euclidean
     * space where bond, bump and separation constraints are not violated.
     * C-alpha coordinates are in Model. In the first overlaid version, 
     * Maxiter is the maximal number of spectral gradient iterations,
     * Eps is the relative stress change.
     * Returns the final stress or a negative float on error.
     * Also sets the Noconv int variable to non-0 if the Spectral Gradient
     * did not converge.
     * In the second overlaid version, the cluster layout is in Pieces
     * and Checkflags tells the routine what to update (uses Willie's
     * simple but efficient "pairwise displacement" optimisation).
     * No value returned.
     */
    float adjust_xyz(Points_& Model, int Maxiter, float Eps, int& Noconv);
    void adjust_xyz(const Trimat_& Dista, Points_& Model,
	const Pieces_& Pieces, int Checkflags) const;
    
	// Auxiliaries
    private:
    
    static float make_iddist(float Actual, float Low, float Up);
    static float limit_iddist(float Ideal, const Restraints_& Restraints, 
	    unsigned int i, unsigned int j);

};
// END OF CLASS Steric_

// ==== GLOBAL PROTOTYPES ====

/* update_bonddist(): re-calculates the first and second
 * squared neighbour distances in Model and puts them into the
 * corresponding off-diagonals of Dista.
 */
void update_bonddist(const Points_& Model, Trimat_& Dista);

// ==== END OF HEADER Steric.h ====

#endif	/* STERIC_CLASSES */
