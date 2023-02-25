#ifndef HOMODEL_CLASS
#define HOMODEL_CLASS

// ==== PROJECT DRAGON: HEADER Homodel.h ====

/* Class for distance-based homology modelling. */

// SGI C++, IRIX 6.2, 16. Jan. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <math.h>

// ---- MODULE HEADERS ----

#include "Polymer.h"
#include "Restr.h"

// ---- UTILITY HEADERS ----

#include "Points.h"
#include "List1.h"
#include "Hirot.h"
#include "Vector.h"
#include "Bits.h"
#include "pdbprot.h"

// ==== CLASSES ====

/* Class Homodel_: implements a distance-based homology modelling
 * approach. Interatomic distance restraints (between CA atoms)
 * are obtained from known structures the sequence of which are
 * aligned to the target sequence and other homologous sequences
 * (see "Align" and "Polymer" for details). The strictness of the
 * restraints is obtained from the pairwise conservation data.
 * Can be asked to provide a restraint list.
 */
class Homodel_
{
    // nested struct
    private:
    
    /* struct Known_: holds information about a known structure
     * which is homologous to the target.
     */
    struct Known_
    {
	Points_ Cas;	// C-alpha coordinates
	char *Seq;	// the sequence
	int Seqidx;	// index of the sequence in the alignment in Pol (-1 if not present)
	float Sim;  // rough similarity to the target sequence
	
	Known_(): Seq(NULL), Seqidx(-1), Sim(0.0) {}
	~Known_() { delete [] Seq; }
    };
    
    // data
    private:
    
    const Polymer_& Pol;  // the associated polymer w/ the alignment and cons data
    Known_ *Knownstructs;   // the array of known structure data
    Known_ *Bestknown;	   // pointer to known struct most similar to target
    unsigned int Knownno;   // no. of known structures
    Bits_ Knownmask, Modelmask;	// equivalence masks for (best) known struct and model
    Vector_ Weight;	// weights for known:model RMS comparison
    Hirot_ Hr;	    // RMS rotation
    
    // methods
    public:
    
	// constructors & destructor
    /* Associates the object with the polymer object P. */
    Homodel_(const Polymer_& P): Pol(P), 
	    Knownstructs(NULL), Bestknown(NULL), Knownno(0) {}
    
    ~Homodel_() { delete [] Knownstructs; }
    
	// access
    // known_no(): returns the number of known structures.
    int known_no() const { return(Knownno); }
    
	// modification
    /* read_knownstr(): reads the known structure(s) from a PDB file Pdbf.
     * More than one structure may be specified in the file as separate
     * chains. The sequences of them are extracted and compared to the
     * sequences already in the alignment: chains which were not found in
     * the alignment are skipped. The C-alpha coordinates are then stored.
     * Nothing is done and a negative integer is returned on error
     * (inaccessible or unreadable PDB file), otherwise the number of
     * structures successfully identified is returned. 0 is returned
     * if NULL or "" was specified as an input PDB file.
     */
    int read_knownstr(const char *Pdbf);
    
    /* NOTE: the following method is used for PVM runs only */
    #ifdef USE_PVM
	/* str_readknown(): reads the PDB structure from the string Pdbstr.
	 * This is a horrible kludge: the string is written to a temporary
	 * file which is then processed by read_knownstr(). The reason is
	 * that the PDB reader is in C, and it wouldn't be easy to convert
	 * all those scanf()-s to sscanf()-s. This method is provided
	 * for the PVM inter-process message passing.
	 * Return values are exactly the same as in read_knownstr().
	 */
	int str_readknown(const char *Pdbstr);
    #endif
    
	// restraint generation
    /* make_restrs(): builds a list of CA distance restraints
     * for all residue pairs in the known structure(s) which participate
     * in the alignment and are closer than the Maxdist (not squared).
     * The lower bound is the shortest distance found
     * in the known structures, the upper bound is the largest.
     * The residue pairs must be separated by Minsepar residues in
     * the sequence (Minsepar==2 is the minimum, allowing for
     * (i, i+2) pairs). 2 is also the default.
     * From V4.11.2 on, also prepares for RMS checks between model
     * structure and the best known structure (to get the right enantiomer).
     */
    List1_<Restr_> make_restrs(float Maxdist, int Minsepar=2);
    
    /* hand_check(): compares the model C-alpha coordinates in Model to
     * the C-alpha coordinates of the scaffold structure most homologous to the model.
     * Returns 1 if Model is more similar to the scaffold than its mirror image, 
     * -1 if a flip was needed (which is done inside) and 0 on error
     * or if Model was not 3D. Cf. the hand_check() fn in "Sterchem".
     */
    int hand_check(Points_& Model);
    
    // hidden methods
    private:

    static void get_ca(const Chain_ *Chain, Points_& Cas);
    
    // forbidden methods
    Homodel_();
    Homodel_(const Homodel_&);
    Homodel_& operator=(const Homodel_&);
    
};
// END OF CLASS Homodel_ 

// ==== END OF HEADER Homodel.h ====
#endif
