#ifndef BETA_CLASS
#define BETA_CLASS

// ==== PROJECT DRAGON: HEADER Beta.h ====

/* Beta-sheet topology and geometry class. Cf. "Segment.h" and
 * "Sstrbase.h" for base class information.
 */

// SGI C++, IRIX 6.2, 14. Aug. 1996. (C) Andras Aszodi

// ---- MODULE HEADERS ----

#include "Segment.h"
#include "Sstrbase.h"

// ==== CLASSES ====

/* Class Beta_ : implements a beta-sheet composed of >=2 strands.
 * Derived from Sheet_ and Sstrbase_ . Can represent intramolecular sheets only
 * and cannot handle bulges, irregularities etc.
 */
class Beta_ : public Sheet_, public Sstrbase_
{
    // data
    protected:
    
    /* Beta-sheets have phasings "up" or "down" and it cannot be known
     * in advance whether a given residue will be "above" or "below" the
     * sheet surface. Ideal sheets with both phasings are generated
     * and stored in Idup and Iddown.
     */
    Points_ Idup, Iddown;
    
    /* The ideal distances will be stored in a matrix: based on "Idup" */
    Trimat_ Dist;
    
    // methods
    public:
    
	// constructors
    /* The ctors call the corresponding Sheet_ ctors and adjust
     * the length of the ideal structure arrays etc.
     */
    Beta_() : Sheet_(), Idup(1, 3), Iddown(1, 3), Dist(1) {}
    Beta_(const Strand_& Str1)
	: Sheet_(Str1), Idup(Str1.len(), 3), 
	    Iddown(Str1.len(), 3), Dist(Str1.len()) {}
    Beta_(const Sheet_& Sh);
    
	// "virtual constructor"
    void operator()(Sstrbase_*&) const;
    
	// dynamic type identification
    bool is_helix() const { return(false); }
    bool is_beta() const { return(true); }
    
	// H-bond topology
    /* hbond_prev(), hbond_next(): return the no. of the previous
     * or the next residue H-bonded to Res or -1 if there's no partner
     * (at sheet edges) or -2 if Res is not a member of the sheet.
     * A warning is also printed in this case.
     */
    int hbond_prev(unsigned int Res) const;
    int hbond_next(unsigned int Res) const;

	// ideal geometry

    /* make_idstruct(): generates two ideal beta-sheets ("up" and "down"
     * phasing) if Changed (inherited from Segmbase_) is true, and puts
     * the 3D coordinates into Idup and Iddown. Returns no. of residues or
     * 0 if something went wrong.
     */
    unsigned int make_idstruct();
    
    /* ideal_dist(): puts the ideal beta-sheet UNsquared distances into
     * a distance matrix Dmat in the right position. Does nothing if
     * the sheet does not fit in the matrix. Prints a warning if Changed==true, 
     * since this indicates that the size was changed without updating the
     * ideal structure and therefore what is returned may be incorrect.
     * (The situation will be dealt with elegantly when the "mutable"
     * keyword finds its way into the compiler.) The ideal distances
     * will be applied at strictness Strict (from Sstrbase_) into Strimat.
     */
    void ideal_dist(Trimat_& Dmat, Trimat_& Strimat) const;
    
    /* ideal_struct(): applies the ideal sheet coordinates stored
     * inside onto the point set Model. Model must be large enough to contain
     * the sheet and when masked with the sheet's mask, the active region
     * must be 3-dimensional. If this is the case, then the ideal structure
     * will be RMS fitted onto Model's active region, the original segment
     * replaced by the rotated/transposed ideal, and the RMS value returned.
     * The phasing will be chosen so that the ideal sheet with a better RMS
     * will be fitted (no knowledge of the phasing is necessary).
     * -1.0 is returned on error. Model's original activation pattern is
     * always retained. Prints a warning if an update is needed.
     */
    double ideal_struct(Points_& Model) const;
    
	// Handedness
    /* check_torsion(): walks over each strand in 3D in Model and calculates all
     * the (i+1, i, k, m) torsion angles where k is i's partner in the next
     * strand, and m is i+1's partner. The angle should be negative.
     * Good and Bad will be set to the no. of correct
     * and incorrect torsion angles.
     * Return value: 1 if Good>=Bad, -1 if Good<Bad, 0 if not in 3D.
     */
    int check_torsion(Points_& Model, 
	    unsigned int& Good, unsigned int& Bad) const;

	// input (for output, see Sstrbase_ )
    /* >>: reads in a sheet description from the input stream In.
     * The format is:-
     * "SHEET [strict]\n"
     * "STRAND <beg> <end>\n"
     * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
     * .....
     * "END\n"
     * The residue position numbers are >=1. The optional [strict]
     * parameter (a floating-point number between 0.0 and 1.0) controls
     * the strictness at which the ideal beta-structure should be 
     * applied. If missing, this strictness is 1.0.
     * Will not update the Beta_ object Beta if errors are detected.
     * The "fail" bit is set on error.
     */
    friend istream& operator>>(istream& In, Beta_& Beta);
    
    // hidden methods
    protected:
    
    void make_ths();
    void write_to(ostream& Out) const;

};
// END OF CLASS Beta_

// ==== END OF HEADER Beta.h ====
#endif	/* BETA_CLASS */
