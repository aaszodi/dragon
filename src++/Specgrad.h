#ifndef SPECGRAD_CLASS
#define SPECGRAD_CLASS

// ==== PROJECT DRAGON: HEADER Specgrad.h ====

/* Majorization algorithm employing the Spectral Gradient Method
 * (Glunt W, Hayden TL, Raydan M, J. Comput. Chem. 14:114-120 (1993)).
 */

// SGI C++ 7.1, IRIX 6.2, 9. May 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <math.h>

// ---- MODULES ----

#include "Points.h"

// ==== CLASSES ====

class Specgrad_
{
    // constants
    static const double SMALL;	// unsafe to divide with things smaller than this
    
    // data
    protected:
    
    Trimat_ Wgt, Distact, Bmat, Smat;   // weight, actual dists, aux matrices
    Matrix_ Xt;	// transposed matrix of coordinates
    Matrix_ Negrad, Oldnegrad; // negative gradients of the stress function
    Matrix_ Xtbest; // coords with the lowest stress
    
    double Wnorm;    // weight matrix norm
    unsigned int N, D;	// no. of active points, dimension
    
    // methods
    public:
    
    /* inits for Size vectors in Dim dimensions */
    Specgrad_(unsigned int Size=10, unsigned int Dim=3):
	Wgt(Size), Distact(Size), Bmat(Size), 
	Smat(Size), Xt(Size, Dim), Xtbest(Size, Dim), 
	Negrad(Size, Dim), Oldnegrad(Size, Dim),
	Wnorm(1.0), N(Size), D(Dim) {}
    
    /* weight(): sets up the calling object to work with a given
     * weight matrix W (with entries >=0.0).
     * Returns the size of the problem.
     */
    int weight(const Trimat_& W);
    
    /* iterate(): performs the iteration on the point set Coords
     * (all vectors are assumed to have the same dimension). The coordinates
     * will be massaged towards the ideal UNsquared distances in Id. 
     * Eps is the relative precision (default 0.001).
     * Itno is the maximal number of iterations which contains
     * the actual iteration number on return. If no convergence was reached
     * then Itno is set to -Itno on return.
     * Return value: the "stress" (weighted dist difference).
     * Negative stress values indicate serious errors.
     */
    float iterate(const Trimat_& Id, Points_& Coords, 
	    int& Itno, float Eps=0.001);
    
    protected:
    
    void make_smat();
    void norm_weights(const Trimat_& Id);
    void actual_dist();
    void make_bmat(const Trimat_& Distid);
    void make_negrad();
    void update_coords(float Alpha);
    void make_alpha(float& Alpha) const;
    float stress(const Trimat_& Distid) const;
};
// END OF CLASS Specgrad_

// ==== END OF HEADER Specgrad.h ====
#endif	/* SPECGRAD_CLASS */
