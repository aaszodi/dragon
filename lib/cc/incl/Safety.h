#ifndef SAFETY_CLASS
#define SAFETY_CLASS

// ==== Safety.h ====

// This tiny class implements a few safety precautions
// which can be switched on and off: safe division that
// traps div-by-zero errors, provides a "small number"
// with which it is unsafe to divide a double, and an improved
// implementation of hypot() for those architectures which
// do not provide it.

// 24-Aug-1998. Andris Aszodi

// SUN standard headers
#if defined(__sun)
#include <limits.h>
#include <float.h>
#endif

/* NOTE: The SGI N32/N64 compilers recognise the built-in bool type
 * and define the macro _BOOL. This workaround is provided for
 * the SGI O32 compiler.
 */
#if !(defined(_BOOL) || defined(HAS_SEEN_BOOL))
#define HAS_SEEN_BOOL
typedef unsigned char bool;
static const bool false=0;
static const bool true=1;
#endif	/* _BOOL */

// ==== CLASSES ====

class Safety_
{
    static const double SMALL;	// smallest non-0 for which 1.0/SMALL is OK
    bool Usesafediv, Usehypot;
    
    // methods
    public:
    
    // ctor: set Safediv to true if you want to use the safe-division
    // feature (default).
    // The Usehypot constant will be initialised in an OS-dependent manner.
    Safety_(bool Safediv=true)
	    : Usesafediv(Safediv), 
	// !!! ---- begin OS-dependent part ---- !!!
	    
	// SGI: hypot() is not strict ANSI, available only
	// if the X/Open compatibility mode is selected
	#if defined(__sgi)
	#if defined(_XOPEN4) && defined(_NO_ANSIMODE)
	    Usehypot(true)
	#else
	    Usehypot(false)
	#endif
	#endif
	
	// Linux/GCC: hypot() available
    	#if defined(__linux__)
	    Usehypot(true)
    	#endif
	
	// SUN
	#if defined(__sun)
	    Usehypot(true)
	#endif
	
	// add any other OS here
    	#if defined(__any_other_OS__)
	    Usehypot(...)
    	#endif
	
	// !!! ---- end of OS-dependent part ---- !!!
	    {}
    
    // return the smallest number which can still divide 1.0 w/o problems 
    double small() const { return SMALL; }
    
    // safe_div(): returns the Usesafediv variable.
    // safe_div(S): sets Usesafediv, returns old value.
    bool safe_div() const { return Usesafediv; }
    bool safe_div(bool S) { bool Olds=Usesafediv; Usesafediv=S; return Olds; }
    
    // no_hypot(): returns true if the Num. Recipes pythag() was to be used
    // no_hypot(H): use Num. Recipes pythag() if H is true, returns old val
    bool no_hypot() const { return !Usehypot; }
    bool no_hypot(bool H) { bool Oldh=!Usehypot; Usehypot=!H; return Oldh; }
    
    // safe_div(): returns Num/Denom if Denom is reliably non-0.
    // If it is closer to 0 then the internal parameter SMALL, then
    // it will be replaced by 100*SMALL and the division is performed
    // with a warning. The line number is printed when greater than 0 (default).
    double safe_div(double Num, double Denom, int Lineno=0) const;

    // pythag(): returns the value sqrt(a^2+b^2) without over-
    // or underflows. On some machines the math library contains
    // the function hypot(x, y) which does the same job and will be used
    // whenever possible. 
    double pythag(double a, double b) const;
};
// END OF CLASS Safety_

// ==== END OF HEADER Safety.h ====

#endif	/* SAFETY_CLASS */
