#ifndef __SPL_CLASS__
#define __SPL_CLASS__

// ====	HEADER Spl.h ====

/* Classic third-order splines. Based on the "spl" module
 * written in C which is based on a Pascal routine collection
 * written in Hungary ages ago...
 */

// SGI C++ 4.0, IRIX 5.3, 10. May 1995. Andris 

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>

// ---- INCLUDE FILES ----

#include "Array.h"

// ---- DEFINITIONS ----

#define SPL_MAX1DER (1e30)

// ==== CLASSES ====

/* Class Spl_ : implements third-order splines. The double arrays
 * X and Y hold the indep. and dep. variable values, Y2 the second
 * derivatives and Yin the values of the integrals from X[0].
 * All arrays have the same length.
 */
class Spl_
{
    // data
    Array_<double> X;	// independent variable
    Array_<double> Y;	// dependent variable
    Array_<double> Y2;	// second derivatives
    Array_<double> Yin;	// integrals
    int Ev;		// 0 if re-evaluation is needed
    
    // methods
    public:
    
	// constructors
    Spl_(unsigned int N=1): X(N), Y(N), Y2(N), Yin(N), Ev(0) {}
    
    // full array access
    const Array_<double>& x_arr() const { return(X); }
    void x_arr(const Array_<double>& Xa);
    const Array_<double>& y_arr() const { return(Y); }
    void y_arr(const Array_<double>& Ya);
    
    // data point access
    double x(unsigned int Idx) const { return(X[Idx]); }
    double& x(unsigned int Idx) { return(X[Idx]); }
    double y(unsigned int Idx) const { return(Y[Idx]); }
    double& y(unsigned int Idx) { return(Y[Idx]); }
    
    // size and data reset
    unsigned int len() const { return(X.len()); }
    void len(unsigned int L) 
    { X.len(L); Y.len(L); Y2.len(L); Yin.len(L); Ev=0; }
    void reset() 
    {
	X.set_values(0.0); Y.set_values(0.0); 
	Y2.set_values(0.0); Yin.set_values(0.0); Ev=0;
    }
    
	// fitting and evaluation
    /* fit_spl(): fits a cubic spline to a series of data points already
     * stored in the calling object. Does not do anything if Ev=1.
     *   yp1,ypn: user-supplied first derivatives at the 1st and last
     *            points. If they are >=SPL_MAX1DER then the second derivatives are
     *            assumed to be 0 (natural spline, default)
     *   Return value: 0 if OK, <0 if something went wrong.
     */
    int fit_spl(double yp1=SPL_MAX1DER, double ypn=SPL_MAX1DER);
    
    /* eval_spl: evaluates a spline at a point xi. If xi is out of
     * range then 0.0 is returned. Otherwise the return value is y(xi). 
     * Adapted from Numerical Recipes.
     * The derivatives and the integral are returned only if
     * they are not NULL: this way,  the user can calc any combination
     * of the first three derivatives and the integral.
     *   xi: interpolation point
     *   Der1..3: 1st..3rd derivative of y at xi
     *   Integ: definite integral of y(x) between X[0] and xi.
     * The derivative and integral ptrs are NULL by default.
     * Prints a warning and returns 0.0 if the values were modified (Ev=0).
     */
    double eval_spl(double xi, double *Der1=NULL, double *Der2=NULL,
		    double *Der3=NULL, double *Integ=NULL) const;

    /* integ_spl(): returns the definite integral of the spline Spl
     * between Low and Up or 0.0 if range error.
     */
    double integ_spl(double Low, double Up) const;

    private:
    /* int_0x(): returns the definite integral of the spline between X[0] and xi.
     * xi is assumed to be within the legal range (this method is called
     * by integ_spl() only and is private).
     */
    double int_0x(double xi) const;
};
// END OF CLASS Spl_

// ==== END OF HEADER Spl.h ====

#endif
