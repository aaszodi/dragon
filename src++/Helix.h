#ifndef HELIX_CLASS
#define HELIX_CLASS

// ==== PROJECT DRAGON: HEADER Helix.h ====

/* Helix topology and geometry class. Cf. "Segment.h" and
 * "Sstrbase.h" for base class information.
 */

// SGI C++, IRIX 6.2, 14. Aug. 1996. (C) Andras Aszodi

// ---- MODULE HEADERS ----

#include "Segment.h"
#include "Sstrbase.h"

// ==== CLASSES ====

/* Class Helix_ : implements a helical segment of a model chain.
 * Derived from Linsegm_ (the linear segment class) and Sstrbase_
 * (the secondary structure abstract base class). 3/10, alpha
 * and pi-helices (all right-handed) are implemented.
 * Helices can be anywhere between positions [1..Rno] on the model chain:
 * positions 0 and Rno+1 are reserved for the N/C-terminal moieties (>=V4.8.1)
 */
class Helix_ : public Linsegm_ , public Sstrbase_
{
    // helix type 
    public:
    enum  Helixtype_ {HX310, ALPHA, HXPI};
    
    // data
    private:
    
    // constants
    static const int HX310_DIAG;  // 3/10-helix H-bond (i,i+2)
    static const int ALPHA_DIAG;  // alpha-helix H-bond (i,i+3)
    static const int HXPI_DIAG;  // pi-helix H-bond (i,i+4)

    static const double RADIUS_310, PITCH_310, TURN_310; // helical params
    static const double RADIUS_ALPHA, PITCH_ALPHA, TURN_ALPHA; 
    static const double RADIUS_PI, PITCH_PI, TURN_PI;
    
    static Array_<double> Dist310;	// ideal distances
    static Array_<double> Distalpha;
    static Array_<double> Distpi;

    Points_ Id;	    // the ideal conformation
    Helixtype_ Htype;	// the helix type
    int Diag;	    // H-bond phase
    
    // methods
    public:
    
	// constructors
    /* Inits the helix to begin at Start and end at Stop to have type Ht
     * (the default is ALPHA).
     * Stop>=Start but the Linsegm_ base class ctor takes care of this.
     * Will be set to store 2 tetrahedral point sets.
     */
    Helix_(unsigned int Start=1, unsigned int Stop=ALPHA_DIAG, Helixtype_ Ht=ALPHA);
	
    /* Inits the helix with a Linsegm_ object with type Ht (default ALPHA). */
    Helix_(const Linsegm_& Ls, Helixtype_ Ht=ALPHA); 
    
	// "virtual constructors"
    void operator()(Sstrbase_*&) const;
    
	// dynamic type identification
    bool is_helix() const { return(true); }
    bool is_beta() const { return(false); }
    int helix_type() const { return(Htype); }
    
	// H-bond topology
    /* hbond_prev(), hbond_next(): return the no. of the previous
     * or the next residue H-bonded to Res or -1 if there's no partner
     * (at helix ends) or -2 if Res is not a member of the helix.
     * A warning is also printed in this case.
     */
    int hbond_prev(unsigned int Res) const;
    int hbond_next(unsigned int Res) const;
    
	// ideal geometry
    /* make_idstruct(): generates the 3D ideal right-handed alpha-helix
     * in Id if the sentinel Changed (inherited from Segmbase_) is true.
     * Returns length or 0 if something fails. 
     */
    unsigned int make_idstruct();
    
    /* ideal_dist(): puts the ideal alpha-helical UNsquared distances into
     * a distance matrix Dmat in the right position. Does nothing if
     * the helix does not fit in the matrix. Prints a warning if Changed==true, 
     * since this indicates that the size was changed without updating the
     * ideal structure and therefore what is returned may be incorrect.
     * (The situation will be dealt with elegantly when the "mutable"
     * keyword finds its way into the compiler.) The adjustment is
     * applied with a variable strictness (Strict is in Sstrbase_ ).
     * Also sets the corresponding strictness entries in Strimat.
     */
    void ideal_dist(Trimat_& Dmat, Trimat_& Strimat) const;
    
    /* ideal_struct(): applies the ideal alpha-helical coordinates stored
     * inside onto the point set Model. Model must be large enough to contain
     * the helix and when masked with the helix's mask, the active region
     * must be 3-dimensional. If this is the case, then the ideal structure
     * will be RMS fitted onto Model's active region, the original segment
     * replaced by the rotated/transposed ideal, and the RMS value returned.
     * -1.0 is returned on error. Model's original activation pattern is
     * always retained. Prints a warning if an update is needed.
     */
    double ideal_struct(Points_& Model) const;
    
	// handedness
    /* check_torsion(): walks over the helix in 3D in Model and calculates all
     * (i, i+3) torsion angles. For right-handed helices these should
     * be all positive. Good and Bad will be set to the no. of correct
     * and incorrect torsion angles.
     * Return value: 1 if Good>=Bad, -1 if Good<Bad, 0 if not in 3D.
     */
    int check_torsion(Points_& Model, 
	    unsigned int& Good, unsigned int& Bad) const;

	// input (for output, see Sstrbase_)
    /* >>: reads in a helix from the stream In. Expects the following
     * format:-
     * "<str> <beg> <end> [strict]\n" 
     * where the spaces represent one or more
     * whitespaces, <beg> and <end> are positive integers. 
     * <str> is one of the following strings: "HX310" for a 3/10 helix, 
     * "ALPHA" or "HELIX" for an alpha-helix, "HXPI" for a pi-helix.
     * [strict] is an optional floating-point number between 0.0 and 1.0, 
     * which, if specified, tells the system to what extent to apply the ideal
     * helical structure in the model.
     * Input stops after <end> or [strict] and the newline is not consumed.
     * On error, the "fail" bit is set and the Helix_ object will not
     * be modified. 
     */
    friend istream& operator>>(istream& In, Helix_& H);
    
    // hidden methods
    protected:
    
    void set_diag();
    void update_iddist(Array_<double>& Dist) const;
    void copy_iddist(Trimat_& Dmat, Trimat_& Strimat, const Array_<double>& Dist) const;
    void make_ths();
    void write_to(ostream& Out) const;

};
// END OF CLASS Helix_

// ==== END OF HEADER Helix.h ====
#endif	/* HELIX_CLASS */
