#ifndef PIECES_CLASS
#define PIECES_CLASS

// ==== PROJECT DRAGON: HEADER Pieces.h ====

/* Class for keeping track of secondary structures and
 * general segments for hierarchic projection.
 */

// SGI C++ 7.1, IRIX 6.2, 4. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>

// ---- UTILITY HEADERS ----

#include "List1.h"
#include "Array.h"
#include "Bits.h"
#include "Smartptr.h"

// ---- MODULE HEADERS ----

/* The following headers contain the declaration of the various
 * secondary structure segment classes in the "pieces" library.
 */
#include "Segment.h"
#include "Sstrbase.h"
#include "Helix.h"
#include "Beta.h"

// ---- TYPEDEFS ----

/* The lists inside Pieces_ hold various derived class objects
 * as smart ptrs (Smartptr_ objects) to their base classes. 
 */
typedef Smartptr_<Sstrbase_> Sstr_ ;

// ==== CLASSES ====

/* Class Pieces_ : holds a list of secondary structures
 * which can be read from a text file and the
 * intervening coil regions are deduced from them.
 * Can be queried for individual sstr elements or segment
 * masks.
 */
class Pieces_
{
    // cluster types
    public:
    enum Clutype_ {UNKNOWN=-1, COIL=0, HELIX, SHEET};	// unknown for invalid clus
    
    // data
    private:
    
    List1_<Sstr_> Secs;	// secondary structure list
    List1_<Linsegm_> Coils;	// coils between sstr elements
    Bits_ Secsmask;	// true for all residues within secondary structures
    Array_<Bits_> Clus;	// mask for each cluster to be used in projection
    Array_<Clutype_> Ctype;	// Ctype[i] is the cluster type of cluster i
    int *Ptclu;    // Ptclu[i] is the cluster index holding i
    unsigned int Rno;	// chain size (may be changed but not queried)
    bool Changed;	// sentinel (true if re-calc needed)
    
    // methods
    public:
    
	// constructors
    /* Inits the object to accept secondary structures within a Resno-long
     * chain. Note that Resno must always be specified explicitly to avoid
     * confusion and therefore there's no default value for it. Also, the
     * default ctor is disallowed (see "forbidden methods").
     */
    Pieces_(unsigned int Resno);
    
	// destructor
    ~Pieces_() { delete [] Ptclu; }
    
	// size
    /* clu_no(): returns the number of clusters. */
    unsigned int clu_no() const { return(Clus.len()); }
    
    /* res_no(): sets the new size to R. Destroys the secstr list and 
     * represents the chain as one long coil (similar to the ctor).
     * The new secstr assignment can then be built by reading a new
     * specification.
     */
    void res_no(unsigned int R);
    
	// structure list access
    const List1_<Sstr_>& secs() const { return(Secs); }
    const List1_<Linsegm_>& coils() const { return(Coils); }
    
	// cluster access
    /* clus(i): returns a const reference to the i-th residue cluster safely. */
    const Bits_& clus(unsigned int i) const { return(Clus(i)); }
    
    /* clu_type(i): returns the type of the i-th residue cluster
     * or UNKNOWN when i is invalid.
     */
    Clutype_ clu_type(int i) const { return((i<0 || i>=Ctype.len())? UNKNOWN: Ctype(i)); }
    
    /* clusters(): returns a const ref to the whole cluster array. */
    const Array_<Bits_>& clusters() const { return(Clus); }
    
    /* member(X): returns the index of the cluster mask of which residue X is 
     * a member or -1 if X is not contained by any of the masks in Clus.
     * members(X, Y): returns the index of the mask which contains both X and Y, 
     * or -1 if X and Y are in different clusters.
     */
    int member(unsigned int X) const { return((X>(Rno+1))? -1: Ptclu[X]); }
    
    int members(unsigned int X, unsigned int Y) const
    {
	register int Xclu=Ptclu[X];
	return((Ptclu[Y]==Xclu)? Xclu: -1);
    }
    
    /* hbond_bits(): returns the Secsmask. */
    const Bits_& hbond_bits() const { return(Secsmask); }

	// input/output
    /* read_secstr(): reads the secondary structure specification from
     * file Secf. If Secf=="" or NULL then the secondary structure
     * layout will be reset to all-coil.
     * Returns 1 if OK, -1 on reset, 0 on error.
     */
    int read_secstr(const char *Secf);

    /* >>: reads the secondary structure specification from a stream In
     * into the Pieces_ object P which is supposed to be pre-initialised.
     * Alpha-helices and beta sheets can be specified according to the
     * following syntax (also cf. the Helix_ and Sheet_ input operator comments):-
     * "<helix> <beg> <end> \n" (where <helix>=["HELIX"|"ALPHA"|"HX310"|"HXPI"] )
     * "SHEET\n"		    (or for beta-sheets)
     * "STRAND <beg> <end>\n"
     * "STRAND <beg> <end> [PAR|ANTI] <this> <other>\n"
     * .....
     * "END\n"
     * Helices and sheets may be freely mixed and comment lines (lines
     * beginning with '#') may be present between the descriptions.
     * Note that comments are not allowed between "SHEET"/"END" pairs.
     * The routine builds a temporary secstr list and updates P only
     * if the input was successful. Checks are made to ensure that
     * all items fit into the chain (based on the internal Rno value)
     * and that there is no overlap between helices and other helices/sheets.
     * Sheet/sheet overlap is allowed but a warning is printed.
     */
    friend istream& operator>>(istream& In, Pieces_& P);
    
    /* <<: lists the secondary structure elements to Out. */
    friend ostream& operator<<(ostream& Out, const Pieces_& P);
    
    // internal consistency maintenance
    private:
    
    void make_coils();
    void make_ptidx();
    
    // "forbidden methods"
    private:
    
    /* Default construction, copying and assignment are forbidden */
    Pieces_();
    Pieces_(const Pieces_&);
    Pieces_& operator=(const Pieces_&);
    
};
// END OF CLASS Pieces_

// ==== END OF HEADER Pieces.h ====

#endif	/* PIECES_CLASS */
