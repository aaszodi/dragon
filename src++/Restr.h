#ifndef RESTR_CLASSES
#define RESTR_CLASSES

// ==== PROJECT DRAGON: HEADER Restr.h ====

/* Stores, smooths and retrieves distance restraints
 * in a uniform way.
 */

// SGI C++, IRIX 6.2, 27. Jan. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

// ---- MODULE HEADERS ----

#include "Sqmat.h"
#include "Trimat.h"
#include "List1.h"
#include "Polymer.h"
#include "Pieces.h"

// Some math libraries lack the float equivalents of the
// usual double functions (eg. float sqrtf(float) in addition
// to double sqrt(double). Define NO_MATHFLOATFUNC in these cases.
#ifdef NO_MATHFLOATFUNC
#define fabsf fabs
#define sqrtf sqrt
#endif

// ==== CLASSES ====

/* Class Restr_: holds the external distance restraints. 
 * Can be asked to return the non-squared or squared restraint
 * between two atoms in two residues with a strictness value.
 */
class Restr_
{
    //data
    private:
    
    String_ Atom1, Atom2;   // the atoms spanning the restraint
    unsigned int Pos1, Pos2;	// the amino acid positions [1..Rno]
    float Low, Up, Low2, Up2, // non-squared and squared lower and upper limits
	Strict;	    // strictness [0..1]
    
    // methods
    public:
    
	// constructor
    Restr_(): 
	Atom1("CA"), Atom2("CA"), Pos1(0), Pos2(0), 
	Low(0.0), Up(0.0), Low2(0.0), Up2(0.0), Strict(1.0) {}
    
	// access
    /* atom(Idx): returns a const ref. to the Idx:th atom string
     * where Idx=[1, 2]. Idx<1 -->1, Idx>2 -->2 w/o warning.
     * atom(Idx, Str): changes the Idx=[1, 2] :th atom string to Str.
     */
    const String_& atom(int Idx) const { return((Idx<=1)? Atom1: Atom2); }
    void atom(int Idx, const String_& Str)
    { if (Idx<=1) Atom1=Str; else Atom2=Str; }
    
    /* pos(Idx): returns the Idx:th position where Idx=[1,2]. As above,
     * Idx<1 -->1, Idx>2 -->2 w/o warning.
     * pos(Idx, Pos): sets the Idx:th restrained position.
     */
    unsigned int pos(int Idx) const { return((Idx<=1)? Pos1: Pos2); }
    void pos(int Idx, unsigned int Pos)
    { if (Idx<=1) Pos1=Pos; else Pos2=Pos; }
    
    /* low(), low2(), up(), up2(): w/o arguments, these methods return
     * the unsquared xxx() or squared xxx2() lower or upper restraints.
     * With a float argument, the restraints can be set simultaneously:
     * a call to low(L) is completely equivalent to a call to low2(L*L).
     * Negative arguments are silently converted to their abs values.
     */
    float low() const { return(Low); }
    void low(float L) { Low=fabsf(L); Low2=Low*Low; }
    float low2() const { return(Low2); }
    void low2(float L2) { Low2=fabsf(L2); Low=sqrtf(Low2); }
    float up() const { return(Up); }
    void up(float U) { Up=fabsf(U); Up2=Up*Up; }
    float up2() const { return(Up2); }
    void up2(float U2) { Up2=fabsf(U2); Up=sqrtf(Up2); }
    
    /* strict(): returns the strictness value.
     * strict(S): sets the strictness value to S. S<0.0 values are
     * converted to 0.0, S>1.0 to 1.0 w/o warning.
     */
    float strict() const { return(Strict); }
    void strict(float S)
    {
	if (S<0.0) S=0.0;
	if (S>1.0) S=1.0;
	Strict=S;
    }
    
	// input,output
    /* >>: reads from the stream In assuming the following format:-
     * "Pos1 Pos2 Lowlim Uplim Strict Atom1 Atom2", where the first five things
     * are the Restr_ data members. Atom1 and Atom2 are atom names corresponding
     * to the PDB conventions,  such as "CA" or "1HB2". Lowercase strings are
     * converted to uppercase but their validity is not checked here. 
     * The string "SCC" may also be specified: it means the side-chain
     * centroid. The maximal string length is 4 characters.
     * Residue numbers are [1..Rno] but the upper limit
     * is not checked here. 
     * 
     * NOTE: In earlier versions (3.x and <=4.7) distance restraints
     * could be specified only between C-alpha and pseudo-C-beta (SCC) atoms.
     * The format was:-
     * "Pos1 Pos2 Lowlim Uplim Strict ABstr", where ABstr corresponded to
     * the 2-char regular expression "[aAbB][aAbB]" with 'a' or 'A' indicating
     * a C-alpha atom, 'b' or 'B' a dummy "C-beta" (SCC). To maintain
     * backward compatibility, this format is also recognised.
     */
    friend istream& operator>>(istream& In, Restr_& R);
    
    /* Output: lists to Out corresponding to the input format. */
    friend ostream& operator<<(ostream& Out, const Restr_& R);
};
// END OF CLASS Restr_

/* Class Restraints_: stores distance restraints and the associated
 * "strictnesses" (weights). Can be asked to perform bounds
 * smoothing and to return (un)squared distance limits.
 */
class Restraints_
{
    // data
    private:
    
    List1_<Restr_> Restrs;  // external distance restraints
    Sqmat_ Lowup, Lowup2;   // unsquared (lower triangle) and squared (upper triangle) restraint limits
    Trimat_ Strict;	    // strictness: 0.0 for non-specific restraints
    Array_<double> Maxsepar;	// maximal residue separation
    unsigned int Size;	    // Rno+2
    
    // public methods
    public:
    
	// constructors
    /* Initialises for Rno-long chains (default 0). */
    Restraints_(unsigned int Rno=0):
	Lowup(Rno+2), Lowup2(Rno+2), Strict(Rno+2), Maxsepar(Rno+2), Size(Rno+2)
    {
	flory_constr();
    }
    
	// access
    /* set_size(): adjusts the size for an Rno-long model chain (that
     * is actually made up of Rno+2 points). Zeroes all internal data
     * if the size has changed.
     * Returns old size.
     */
    unsigned int set_size(unsigned int Rno);
    
    // Const access is provided to the external restraint list
    const List1_<Restr_>& ext_restr() const { return(Restrs); }
    unsigned int restr_no() const { return(Restrs.len()); }
    
    // Const access to the maximal allowable separation
    double max_separ(unsigned int S) const
    {
	return(Maxsepar[(S>=Size? Size-1: S)]);
    }
    
    /* low(),low2(),up(),up2(): retrieve or set the unsquared or squared
     * lower or upper limits between the ith and jth C-alphas.
     * Calling low(i, j, L) is equivalent to a call
     * to low2(i, j, L*L). In general, i, j is in the range [1..Rno]
     * and j may be larger or smaller than i. i==j returns 0.0.
     * Negative distance limits are silently converted to their abs values.
     * On index range error, no "set" operations are done and the "get"
     * methods return 0.0 with a warning.
     */
    double low(unsigned int i, unsigned int j) const;
    void low(unsigned int i, unsigned int j, double Low);
    double low2(unsigned int i, unsigned int j) const;
    void low2(unsigned int i, unsigned int j, double Low2);
    double up(unsigned int i, unsigned int j) const;
    void up(unsigned int i, unsigned int j, double Up);
    double up2(unsigned int i, unsigned int j) const;
    void up2(unsigned int i, unsigned int j, double Up2);

    /* strict(): gets or sets the strictness value for the i:j C-alpha restraint.
     * Negative strictnesses are silently converted to absolute values.
     * Returns the correct values for the non-specific restraints
     * (bonds and bumps).
     */
    double strict(unsigned int i, unsigned int j) const;
    void strict(unsigned int i, unsigned int j, double Str)
    {
	if (check_index(i, j)) Strict[i][j]=fabs(Str);
    }
    
    /* specific(): returns "true" if the i,j:th restraint is specific,
     * ie. either comes from a secondary structure assignment or external
     * restraint list. Returns "false" for unspecific restraints such as
     * bond lengths, bumps and maximal separations.
     */
    bool specific(unsigned int i, unsigned int j) const
    {
	return(bool(check_index(i, j) && Strict[i][j]>0.0));
    }
    
    /* hard(): a restraint is hard either if it is specific or if it
     * is a CA:CA bond or geminal distance (1st or 2nd diagonal). 
     */
    bool hard(unsigned int i, unsigned int j) const
    {
	return(bool(abs(i-j)<=2 || specific(i, j)));
    }
    
	// restraint setup
    /* setup_restr(): clears the restraint matrices and sets them up
     * according to the list of external restraints (which should already
     * be prepared in the calling object), secondary structure (from Pieces)
     * and intra-monomer atom distances (from Polymer). Call only once
     * before the simulations.
     */
    void setup_restr(const Pieces_& Pieces, const Polymer_& Polymer);
    
    /* init_distmat(): produces a random (squared) distance matrix with entries from a
     * Gaussian distribution. The average and S.D is calculated from the
     * expected radius of the molecule, assuming a spherical shape.
     * The random number generator is initialised with
     * Randseed. Only random numbers falling between the distance bounds
     * are accepted. "Softly" restrained positions are modified by the
     * hydrophobic distance estimates.
     */
    void init_distmat(Trimat_& Dist, Polymer_& Polymer, 
	    long Randseed) const;
    
    /* exp_rad(): returns the expected radius for an Rno-long
     * chain with a density Dens.
     */
    static float exp_rad(unsigned int Rno, float Dens=0.00636)
    {
	return(pow((3.0F*Rno/(4.0F*M_PI*Dens)), 0.333333));
    }
    
	// input-output
    /* Output: lists to Out the restraints. */
    friend ostream& operator<<(ostream& Out, const Restraints_& D);
    
    /* read_restrs(): reads restraints from a file Fname. Assumes that
     * the file contains comments (lines beginning w/ '#') and
     * lines corresponding to the restraint format (see above in the
     * Restr_ class input for details), one line per restraint each.
     * If Fname==NULL or "" then the restraint list is cleared.
     * Checks the restraints for range conformance (residue nos in the
     * range [1..Rno]) (Rno comes from the length of the Maxsepar array internally.)
     * If Rno==0, then nothing is done (no restraints would pass).
     * Interatom restraints will be converted into CA and SCC
     * (sidechain centroid)-based restraints using the polymer 
     * description in Polymer and the atom:(CA|SCC) distances
     * in Acdist.
     * The old restraint list will be replaced with the new.
     * Return value: 0 on error, 1 if OK.
     */
    int read_restrs(const char *Fname, const Polymer_& Polymer);
    
    /* Input: reads restraints (and comment lines beginning with '#')
     * from the input stream Inf. See >>Restr_ for the restraint line
     * format. The old restraint list is cleared. The new list after
     * input will contain a mixture of interatomic and CA:SCC restraints.
     * To convert to an all-CA:SCC list, invoke convert_restraints()
     * after >> (this is done automagically in read_restrs() above).
     * Note that the density is not set here!
     */
    friend istream& operator>>(istream& Inf, Restraints_& Restraints);
    
    /* add_restrs(): adds the CA|SCC restraints specified in the list
     * Rs to the internal restraint list. Filters out those restraints
     * which are between side-chain atoms. 
     * Useful for restraints coming from the homology
     * modelling module. Call after read_restrs() since that method will
     * clear the restraint list.
     * Return value: the number of restraints processed.
     */
    int add_restrs(const List1_<Restr_>& Rs);
    
    /* convert_restraints(): converts the restraints in the internal restraint list
     * which may have been specified between side-chain atoms, into
     * (looser) restraints between C-alpha atoms (CA) and/or side
     * chain centroids (SCC). Currently DRAGON uses a reduced
     * chain representation with CAs and SCCs only.
     */
    void convert_restraints(const Polymer_& Polymer);
    
    // private methods
    private:
    
    int check_index(unsigned int& i, unsigned int& j) const;
    int merge_restr(unsigned int i, unsigned int j, 
	    double Low, double Up, double Str);
    void setup_bondbump();
    void setup_extrestr(const Polymer_& Polymer);
    void setup_secstrestr(const Pieces_& Pieces);
    int smooth_restr(unsigned int Pass=1);
    void flory_constr();
    static int get_cascc(const Polymer_& Polymer, unsigned int Pos, 
	const String_& Atom, float& Cad, float& Sccd);
    static int absorb_restraint(List1_<Restr_>& Rlist, const Restr_& R);
    void add_restraint(const Restr_& R);
    
    // static members
    public:
    
    static const float CA_BONDANGLE;   // max. CA:CA:CA virtual bond angle (radians)
    static const float CA_BUMP;	// C-alpha bump radius
    static const float CA_1_MIN;	// CA:CA virtual bond length
    static const float CA_1_MAX;
    static const float CA_2_MIN;	// second neighbour distance
    static const float CA_2_MAX;
    
    static const float NTCA_BUMP; // terminal N bump+CA_BUMP
    static const float NT_BONDLEN;  // H3N-C bond length
    static const float NT_2_MIN;    // min. Nt:second CA dist
    static const float NT_2_MAX;    // max. Nt:second CA dist
    static const float CTCA_BUMP; // terminal C bump+CA_BUMP
    static const float CT_BONDLEN;  // C-COO bond length
    static const float CT_2_MIN;    // min. Ct:second CA
    static const float CT_2_MAX;    // max. Ct:second CA
    static const float NTCT_BUMP; // terminal N and C bump together

    static const float STR1;	// CA:CA neighbour strictness
    static const float STR2;	// CA:CA second neighbour strictness
    static const float STRA;	// general CA:CA vdW bump strictness
    static const float STRB;	// strictness for bumps involving sidechains
};
// END OF CLASS Restraints_

#endif	/* RESTR_CLASSES */
