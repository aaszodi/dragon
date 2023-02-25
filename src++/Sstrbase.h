#ifndef SSTRBASE_CLASS
#define SSTRBASE_CLASS

// ==== PROJECT DRAGON: HEADER Sstrbase.h ====

/* Base class for handling secondary structures: H-bond topology
 * and ideal geometry. Chain topology comes from the Segment module.
 */

// SGI C++, IRIX 6.2, 14. Aug. 1996. (C) Andras Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "Bits.h"
#include "Points.h"
#include "Array.h"
#include "Trimat.h"

// ---- MODULE HEADERS ----

#include "Segment.h"

// ==== CLASSES ====

/* struct Thidx_: stores the indices of 4 points which define a
 * tetrahedron for secstr de-tangling (cf. "Tangles" module).
 * Helices will have in general 2 tetrahedra superimposed on them, 
 * sheets will have one less than the number of strands.
 */
struct Thidx_
{
    unsigned int P1, P2, P3, P4;    // tetrahedron point indices
    
    // default ctor: demanded by Array_<Thidx_>
    Thidx_(): P1(0), P2(0), P3(0), P4(0) {}
};

/* Class Sstrbase_ : an abstract base class that serves as the ancestor
 * of the Helix_ and Beta_ classes. Contains pure virtual fns 
 * declaring the methods for obtaining H-bond partnership information
 * and ideal geometry data. Can be asked for making an ideal helix
 * (needed for beta-strands as well).
 */
class Sstrbase_ : public virtual Segmbase_
{
    // data
    protected:
    Array_<Thidx_> Thedra;  // tetrahedron indices for detangling
    float Strict;   // strictness value for ideal structure enforcement
    
    // methods
    public:
    
	// constructor
    /* Inits to hold Thno tetrahedra: default is 0 */
    Sstrbase_(unsigned int Thno=0): Thedra(Thno), Strict(1.0) {}
    
	// "virtual constructor"
    /* Overloads the fn call operator to implement a "virtual ctor"
     * which creates an instance of the derived objects and returns
     * a base class ptr to them. Needed by the Smartptr_ template class.
     */
    virtual void operator()(Sstrbase_*& Sptr) const =0;
    
	// dynamic type identification
    /* is_helix(), is_beta(): these methods go against the orthodox
     * C++ polymorphism philosophy: however, I need them for the
     * conversion to C structs (in Output.c++).
     */
    virtual bool is_helix() const =0;
    virtual bool is_beta() const =0;
    
	// H-bond topology
    /* hbond_prev(), hbond_next(): return the no. of the previous
     * or the next residue H-bonded to Res or -1 if there's no partner
     * (at helix ends or sheet edges) or -2 if Res is not a member of
     * the structure.
     */
    virtual int hbond_prev(unsigned int Res) const =0;
    virtual int hbond_next(unsigned int Res) const =0;
    
	// Tetrahedral points
    /* get_thedra(): returns a const ref. to the tetrahedron index array. */
    const Array_<Thidx_>& get_thedra() const { return(Thedra); }
    
	// ideal geometry
    /* NOTE: updating the internals of the ideal geometry is
     * time-consuming. It should be done whenever the internal sentinel
     * Changed==true. The ideal_dist2() and ideal_struct() methods
     * should check for this condition and invoke make_idstruct()
     * when needed, thus making the update transparent. However, 
     * in the current environment they cannot be declared then "const".
     * The solution would be to declare Changed, the ideal structure
     * internals and the related methods "mutable" so that the
     * outside "const"-ness is retained and yet updating is deferred
     * to only when it is needed. I cannot do this with the present
     * compiler so the workaround is to update the secondary structure
     * objects manually in "Pieces.c++" after input. 
     */
    
    /* make_idstruct(): generates the 3D ideal structure in the calling
     * object if the sentinel Changed (inherited from Segmbase_) is true.
     * Returns length or 0 if something fails. 
     */
    virtual unsigned int make_idstruct() =0;
    
    /* ideal_dist2(): puts the ideal UNsquared distances into
     * a distance matrix Dmat in the right position and the
     * associated strictness in Strict. Does nothing if
     * the structure does not fit in the matrix. 
     */
    virtual void ideal_dist(Trimat_& Dmat, Trimat_& Strict) const =0;

    /* ideal_struct(): applies the ideal secondary structure coordinates stored
     * inside onto the point set Model. Model must be large enough to contain
     * the structure and when masked with the mask, the active region
     * must be 3-dimensional. If this is the case, then the ideal structure
     * will be RMS fitted onto Model's active region, the original segment
     * replaced by the rotated/transposed ideal, and the RMS value returned.
     * -1.0 is returned on error. Model's original activation pattern is
     * always retained.
     */
    virtual double ideal_struct(Points_& Model) const =0;
    
	// Handedness
    /* check_torsion(): checks suitably chosen torsion angles in Model
     * corresponding to the type of secstr. Tests handedness by
     * counting Good and Bad torsion angles.
     * Return value: 1 if Good>=Bad, -1 if Good<Bad, 0 if not in 3D.
     */
    virtual int check_torsion(Points_& Model, 
	    unsigned int& Good, unsigned int& Bad) const =0;
    
	// output
    friend ostream& operator<<(ostream& Out, const Sstrbase_& S);
    
    // hidden methods
    protected:

    virtual void make_ths() =0;
    static unsigned int make_helix(Points_& Hel, double Radius, 
	    double Pitch, double Turn, int Phasing=1);
    static double pos4_angle(const Vector_& P1, const Vector_& P2,
	    const Vector_& P3, const Vector_& P4);
    virtual void write_to(ostream& Out) const =0;

};
// END OF CLASS Sstrbase_

// ==== END OF HEADER Sstrbase.h ====
#endif	/* SSTRBASE_CLASS */
