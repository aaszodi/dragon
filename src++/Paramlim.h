#ifndef PARAMLIM_TMPL_DECLS
#define PARAMLIM_TMPL_DECLS

// ==== PROJECT DRAGON: TEMPLATE HEADER Paramlim.h ====

/* Stores global numeric parameters. */

// SGI C++ 7.1, IRIX 6.2, 7. Mar. 1997. Andris Aszodi

// ---- UTILITY HEADERS ----

#include "Parambase.h"

// ==== CLASSES ====

/* Class Paramlim_ : stores parameters which must fall between a lower
 * and upper bound. Intended for numerical parameters. Derived from Parambase_
 */
template <class T_>
class Paramlim_ : public Parambase_
{
    // data
    protected:
    
    T_ Value, Default;	// the value and its default
    T_ Low, Up;	    // lower and upper limits
    
    // methods
    public:
    
	// constructors
    /* Default ctor: it is assumed that type T_ has a default ctor. */
    Paramlim_(): Parambase_() {}
    
    /* Inits to hold the default value Defval which falls between
     * the lower and upper bounds L and U. The name and description
     * strings are specified in Nm and Ds (both defaults to NULL).
     * If L>U then the values are silently swapped. If Defval is 
     * outside the range defined by [L..U], then it is appropriately
     * modified w/o warning. 
     */
    Paramlim_(const T_& Defval, const T_& L, const T_& U, 
	    const char *Nm=NULL, const char *Ds=NULL);

	// access
    operator T_ () const { return(Value); }
    void reset_default() { Value=Default; Changed=true; }
    
    /* set_deflims(): resets the default value and the limits
     * (similar to the ctor). The value will be reset to the new
     * default value.
     */
    void set_deflims(const T_& Defval, const T_& L, const T_& U);
    
	// input/output
    /* read_from(): reads the following line from an input stream In:
     * "NAME value"
     * where NAME should match the Name member of the object and value
     * should conform to the actual type of Value (T_). If "value" is
     * outside the limits set by Low and Up and it will be modified
     * silently. 
     */
    int read_from(istream& In);
    
    /* write_to(): writes the calling object to the output stream Out
     * in the following format:-
     * "# DESCRIPTION_STRING (default DEFAULT, limits LOW..UP)\n"
     * "NAME value\n"
     * which can be parsed by the input routine. If Comments=false, 
     * the # line is omitted.
     */
    ostream& write_to(ostream& Out, bool Comments=true) const;
    
};
// END OF CLASS Paramlim_

#ifdef INCLUDE_TMPL_DEFS
#include "Paramlim.c++"
#endif

// ==== END OF TEMPLATE HEADER Paramlim.h ====

#endif	/* PARAMLIM_TMPL_DECLS */
