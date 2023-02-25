#ifndef TANGLES_CLASS
#define TANGLES_CLASS

// ==== PROJECT DRAGON: HEADER Tangles.h ====

/* Secondary-structure-based tangle detection and elimination. */

// SGI C++ 4.0, IRIX 5.3, 27. Sept. 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

// ---- MODULE HEADERS ----

#include "Pieces.h"

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "Bits.h"
#include "List1.h"
#include "Points.h"
#include "Svd.h"

// ==== CLASSES ====

/* Class Tangles_: performs secondary-structure-based tetrahedral
 * tangle detection and adjustment. 
 */
class Tangles_
{
    // data
    private:
    
    // struct for storing cluster index pairs which are entangled
    struct Violpair_ 
    {
	unsigned int Idx1, Idx2;
	
	// ctor demanded by Array_<Violpair_>
	Violpair_(): Idx1(0), Idx2(0) {}
    };
    
    List1_<Violpair_> Viols;   // stores indices of entangled segments
    Svd_ Thsvd;	// SVD decomposition of tetrahedra
    Points_ Displ, Ctrs;	    // displacement and centroid vectors
    Array_<unsigned int> Dnos;	// number of displacements
    Bits_ Tmask;		// true for entangled segments (keep centroids)

    // methods
    public:
    
	// constructors
    /* Inits so that tangle checks can be performed on the segments stored
     * in Pieces. (Note that the default ctor is declared private and undefined.)
     */
    Tangles_(const Pieces_& Pieces):
	Viols(), Thsvd(3, 3), 
	Displ(Pieces.clu_no()), Ctrs(Pieces.clu_no()),
	Dnos(Pieces.clu_no()), Tmask(Pieces.clu_no()) {}
    
	// update
    /* update_pieces(): updates the internal data members so that tangle
     * detection can be carried out on the segments in Pieces now.
     * Must be called whenever Pieces changes (no automatic consistency yet).
     */
    void update_pieces(const Pieces_& Pieces);
	
	// tangle detection and elimination

    /* tangle_detect(): checks whether the structure in Xyz is entangled,
     * given the segment layout in Pieces. (The calling object is assumed
     * to have called the update_pieces() method after the last change
     * to Pieces.)
     * Return value: 1 if there was(were) tangle(s), 0 if OK.
     */
    unsigned int tangle_detect(const Pieces_& Pieces, Points_& Xyz);

    /* tangle_elim(): checks and optionally adjusts tangled conformations.
     * The check is carried out by testing whether parts of the chain
     * penetrate tetrahedra put on secondary structure segments.
     * The adjustment is carried out Iter times which
     * moves the entangled segments away from each other (scaled by Tadj).
     * Pieces holds the secondary structure and intervening coil list, 
     * Xyz holds the coordinates and its masking status will be restored
     * upon return.
     * Return value: the no. of entanglements (0 means OK). Also the
     * number of iterations actually done is returned in Iter.
     */
    unsigned int tangle_elim(const Pieces_& Pieces, Points_& Xyz,
	    double Tadj, unsigned int& Iter);

    // hidden methods
    private:
    
	// "forbidden" ctor
    Tangles_();	    // undefined
    
	// tangle detection and elimination

    unsigned int find_tangles(const Pieces_& Pieces, Points_& Xyz, int Adjust=0);
    void adjust_tangles(const Pieces_& Pieces, Points_& Xyz, double Tadj);
    int contain_segment(const Bits_& Segmask, const Points_& Xyz, unsigned int Oidx) const;

	// tetrahedra
    int make_thedron(const Points_& Xyz, const Thidx_& Thidx);
    void make_svect(const Vector_& Vec, const Vector_& Orig, Vector_& S) const;
    int th_viol(const Vector_& Sprev, const Vector_& Snext) const;
    
};
// END OF CLASS Tangles_

// ==== END OF HEADER Tangles.h ====
#endif	/* TANGLES_CLASS */
