// ==== METHODS Safety.c++ ====

// This tiny class implements a few safety precautions
// which can be switched on and off: safe division that
// traps div-by-zero errors, provides a "small number"
// with which it is unsafe to divide a double, and an improved
// implementation of hypot() for those architectures which
// do not provide it.
/* ---- HISTORY ----
 *NWD - 2/7/98 added float.h to includes 
 */
// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <float.h> 
// ---- MODULE HEADER ----

#include "Safety.h"

// ==== Static initialisation ====

const double Safety_::SMALL=sqrt(DBL_MIN)/DBL_EPSILON;	// unsafe to do 1.0/SMALL

// ==== Safety_ MEMBER FUNCTIONS ====

// safe_div(): returns Num/Denom if Denom is reliably non-0.
// If it is closer to 0 then the internal parameter SMALL, then
// it will be replaced by 100*SMALL and the division is performed
// with a warning. The line number is printed when greater than 0 (default).
double Safety_::safe_div(double Num, double Denom, int Lineno) const
{
    if (Usesafediv && fabs(Denom)<SMALL)  // problem
    {
	cerr<<"\n! Svd_::safe_div("<<Num<<", "<<Denom<<"): Dangerous division";
	if (Lineno) cerr<<" at line "<<Lineno<<endl;
	else cerr<<endl;
	Denom=Denom>=0.0? 100*SMALL: -100*SMALL;
    }
    return(Num/Denom);
}
// END of safe_div()

// pythag(): returns the value sqrt(a^2+b^2) without over-
// or underflows. On some machines the math library contains
// the function hypot(x, y) which does the same job and will be used
// whenever possible.
double Safety_::pythag(double a, double b) const
{
    if (Usehypot) return(hypot(a, b));
    else
    {
	register double at=fabs(a), bt=fabs(b), ct;

	if (at<SMALL) return(bt);
	if (bt<SMALL) return(at);
	if (at>bt)
	{
	    ct=bt/at; return(at*sqrt(1.0+ct*ct));
	}
	else
	{
	    ct=at/bt; return(bt*sqrt(1.0+ct*ct));
	}
    }
}
/* END of pythag() */

// ==== END OF METHODS Safety.c++ ====
