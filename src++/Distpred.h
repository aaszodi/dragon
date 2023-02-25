#ifndef DISTPRED_CLASS
#define DISTPRED_CLASS

// ==== PROJECT DRAGON: HEADER Distpred.h ====

/* Interresidue distance prediction based on the conserved
 * hydrophobicity score.
 */

// SGI C++ 4.0, IRIX 5.3, 12. May 1995. Andras Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- UTILITY HEADERS ----

#include "Vector.h"
#include "Spl.h"
#include "Cdf.h"
#include "Array.h"

// ==== CLASSES ====

/* Class Distpred_ : stores the ideal C-alpha distance distribution
 * in a Spl_ (spline) object and can generate a transform function
 * from the conserved hydrophobicity scores in the sequence.
 * Can be queried for an estimated interresidue distance.
 */
class Distpred_
{
    // data
    static Spl_ Idspl;	    // ideal C-alpha CDF approximated by a spline
    Vector_ Par;	    // transform function parameters
    
    // methods
    public:
    
	// constructor
    /* Inits the parameter vector to the default values. */
    Distpred_() { Par=init_par(); }
    
	// initialisation
    /* init_idspl(): fills up the ideal C-alpha distance CDF spline with
     * the observed C-alpha distance distribution of a bunch of monomeric
     * proteins between 100 and 200 residues. See Aszodi & Taylor, J. Math. Chem.
     * for details.
     */
    static Spl_ init_idspl();
    
    /* init_par(): inits the parameter vector to the parameter
     * values in the J. Math. Chem. paper.
     */
    static Vector_ init_par();
    
	// transform function parameter estimation
    /* estim_params(): estimates the parameters of the hydrophobic score
     * transform function from the phobicity*conservation data in
     * Consphob. If Consphob is empty then no action is taken.
     */
    void estim_params(const Array_<double>& Consphob);

    /* dist_phob(H,P): returns the distance corresponding to the
     * raw hydrophobic estimate H or 0.0 if H<=0.0. It is assumed
     * that all P[i]>=0.0.
     * dist_phob(H): the non-static version simply calls the static
     * version with Par as the parameter vector.
     */
    static double dist_phob(double H, const Vector_& P);
    double dist_phob(double H) { return(dist_phob(H, Par)); }
    
    // private methods
    private:
    
    static Cdf_ make_distr(const Array_<double>& Consphob);
    static double transform_hdist(double H, const Vector_& P);
    
    // "forbidden methods"
    Distpred_(const Distpred_&);
    Distpred_& operator=(const Distpred_&);
    
};
// END OF CLASS Distpred_

// ==== END OF HEADER Distpred.h ====
#endif /* DISTPRED_CLASS */
