#ifndef ACCESS_CLASS
#define ACCESS_CLASS

// ==== PROJECT DRAGON: HEADER Access.h ====

/* Conic accessibility calculations. */

// SGI C++, IRIX 6.2, 25. Nov. 1996. Andris Aszodi

// ---- STANDARD HEADERS ---- 

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- MODULE HEADERS ----

#include "Fakebeta.h"
#include "Polymer.h"

// ---- UTILITY HEADERS ----

#include "Bits.h"
#include "Trimat.h"
#include "Points.h"

// ==== CLASSES ====

/* Class Access_ : stores the relative "cone" accessibility of a
 * structure. Can be asked to adjust accessibilities in Eucl space.
 */
class Access_
{
    // data
    private:
    
    /* The shieldedness states are as follows:-
     *	VERY_EXPOSED: <10 %
     *	MED_EXPOSED: 10..15 %
     *	SLGT_EXPOSED: 15..20 %
     *	AVERAGE: 20..80 %
     *	SLGT_BURIED: 80..85 %
     *	MED_BURIED: 85..90 %
     *	VERY_BURIED: >90 %
     */
    enum Shstate_ {VERY_EXPOSED=0, MED_EXPOSED, SLGT_EXPOSED, AVERAGE, 
		SLGT_BURIED, MED_BURIED, VERY_BURIED};
    
    Array_<double> Relsh;   // the relative shield
    
    Array_<double> Di0, Dik;	// centroid dist sq and i-k dist sq
    Array_<unsigned int> Close;	// index lookup
    
    Bits_ Surface, Buried;	    // residues with known accessibilities
    
    // methods
    public:
    
	// constructors
    /* Inits to keep track of Rno>=0 (default 0) amino acids. */
    Access_(unsigned int Rno=0): 
	Relsh(Rno), Di0(2*(Rno+2)), Dik(2*(Rno+2)), 
	Close(2*(Rno+2)), Surface(Rno), Buried(Rno) {}

	// size
    /* set_size(): resets the size of the internal arrays to Rno,
     * the number of residues in the chain. Returns old size.
     * This routine must be called when the chain size changes.
     */
    unsigned int set_size(unsigned int Rno);
	
	// adjustments and scoring 
    /* solvent_xyz(): calculates the accessibility of the structure Xyz
     * given in Euclidean coordinates and performs an accessibility
     * adjustment. The sequence is taken from Polymer. Non-H-bonded
     * residues are moved outwards from the core so they need special 
     * treatment: Hbond is a [0..Rno-1] bit-vector in which every residue participating
     * in a H-bond is marked true (can be obtained from a Pieces_ object).
     * Dim mismatches are checked.
     * Return value: 0 on error, 1 if OK.
     */
    int solvent_xyz(const Polymer_& Polymer, const Bits_& Hbond, 
	    Points_& Xyz);
    
    /* score_dist(),score_xyz(): calculate an accessibility score
     * (optimum 0.0) for distance spaces and Euclidean objects, 
     * respectively.
     * Return a very high value on error.
     */
    float score_dist(const Polymer_& Polymer, const Trimat_& Dista);
    float score_xyz(const Polymer_& Polymer, const Points_& Xyz);
    
	// Input/output
    /* read_file(): reads known accessibilities from a file called Accfnm.
     * If it is NULL (the default) or "", then the surface/buried bit-vectors
     * are cleared. If Accfnm cannot be opened, then nothing is changed.
     * The file format is as follows:-
     * Lines starting with '#' are comments, 
     * [sS] <int> <int> ....\n are the nos. of residues known to be on the surface, 
     * [bB] <int> <int> ....\n are the nos. of residues known to be buried.
     * Any number of residues may be listed after the initial s, S, b, B character
     * provided the total length of the line does not exceed 132 characters,  
     * and any number of lines are allowed. The [sSbB] character must be the
     * first character on the line and the numbers are separated by whitespaces.
     * Duplicates are ignored. If a residue
     * is specified as surface and buried at the same time, it will be regarded
     * as an ordinary residue. Out-of-range residue numbers (outside 1..Rno)
     * are ignored with a warning.
     * Return value: number of lines processed.
     */
    unsigned int read_file(const char *Accfnm=NULL);
    
    /* >>: Reads the list of surface/buried residues from the stream Inf.
     * See also the comments for read_file().
     */
    friend istream& operator>>(istream& Inf, Access_& Acc);
    
    /* <<: lists the residues marked as surface/buried to Out nicely. */
    friend ostream& operator<<(ostream& Out, const Access_& Access);
    
    private:
    void betacone_shield(const Trimat_& Dista, const Polymer_& Polymer);
    int betacone_xyz(const Polymer_& Polymer, const Points_& Xyz);
    float get_score(const Polymer_& Polymer) const;
    static int get_shield(double Rsh, char Aa);
    
};
// END OF CLASS Access_

// ==== END OF HEADER Access.h ====
#endif /* ACCESS_CLASS */
