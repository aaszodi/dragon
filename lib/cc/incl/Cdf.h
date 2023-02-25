#ifndef __CDF_CLASS__
#define __CDF_CLASS__

// ==== HEADER Cdf.h ====

/* Estimation of cumulative distribution functions. */

// SGI C++ 4.0, IRIX 5.3, 11. May 1995. Andris Aszodi 

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "Vector.h"

// ==== CLASSES ====

/* Class Cdf_: the class can be filled up with a set of doubles
 * and the cumulative distribution function can be approximated
 * as a histogram. 
 */
class Cdf_
{
    // data
    private:
    Array_<double> X;	// indep. variable in ascending order
    Array_<double> Y;	// the CDF distribution (all arrays have the same length)
    Array_<unsigned int> Yint;	// the frequency count
    unsigned int N;    // N: total no. of values added
    char Ev;	// evaluation indicator: 0 if data are "dirty"
    
    // methods
    public:
    
	// constructor
    /* Inits to have Binno>=2 bins equally spaced between Low and Up.
     * If Low>Up then they're swapped silently. If Binno<=1 then Binno==2 is used.
     */
    Cdf_(unsigned int Binno, double Low, double Up);
	
	// reset 
    /* reset(): resets the calling object to have Binno>=2 bins equally spaced
     * between Low and Up. Zeroes everything. Very similar to what the ctor does.
     */
    void reset(unsigned int Binno, double Low, double Up);
	
	// data update
    /* Cdf+=V, Cdf-=V: add or remove a value V to/from the Cdf object. */
    Cdf_& operator+=(double V);
    Cdf_& operator-=(double V);
	
	// access
    unsigned int bin_no() const { return(X.len()); }
    const Array_<double>& x_arr() const { return(X); }
    const Array_<double>& y_arr();
    Vector_ x_vec() const;
    Vector_ y_vec();
	
	// private methods
    unsigned int get_index(double V) const;
    int eval_cdf();
	
};
// END OF CLASS Cdf_ 

// ==== END OF HEADER Cdf.h ====
#endif
