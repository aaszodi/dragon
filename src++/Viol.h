#ifndef VIOL_CLASSES
#define VIOL_CLASSES

// ==== PROJECT DRAGON: HEADER Viol.h ====

/* For keeping track of restraint violations. */

// SGI C++ 7.1, IRIX 6.2, 2. Apr. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <math.h>

// ---- MODULES ----

#include "String.h"
#include "List1.h"

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this. Get around by the following macro
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#define fabsf fabs
#endif

// ==== CLASSES ====

/* Viol_: stores the data of the violating residue pair,
 * the type of violation and the extent.
 */
class Viol_
{
    /* type of violations: BOND is for covalent bonds (CA:CA or CA:SCC),
     * NONBD is for unbonded and unrestrained atom pairs, RESTR for
     * externally restrained atom pairs, HELIX and SHEET for 2ostr elements, 
     * UNDEF for anything else
     */
    public:
    enum Violtype_ {UNDEF=0, BOND, NONBD, RESTR, HELIX, SHEET};

    // data
    protected:
    
    Violtype_ Vtype;	    // type of violation
    int Res1, Res2;	    // residue indices [1..Rno]
    String_ Atom1, Atom2;   // "CA" or "SCC"
    float Viol, Strict, Ideal, Actual;    // violation, strictness, ideal and actual distances

    // methods
    public:
	// ctor
    Viol_():
	    Vtype(UNDEF), Res1(0), Atom1("CA"), Res2(0), Atom2("CA"), 
	    Viol(0.0F), Strict(1.0F), Ideal(0.0F), Actual(0.0F) {}
    
    // set the violation type
    void viol_type(Violtype_ V) { Vtype=V; }
    
    // set the atom ID and residue no: Onetwo<=1 for Atom1, Onetwo>=2 for Atom2
    void atom(int Onetwo, const char *Atomnm, int Resno, Violtype_ V=UNDEF)
    {
	if (Onetwo<=1) { Atom1=String_(Atomnm); Res1=Resno; }
	else { Atom2=String_(Atomnm); Res2=Resno; }
	if (V!=UNDEF) Vtype=V;
    }
    
    /* rel_viol(): calculates and stores the weighted relative distance violation.
     * Act is the actual (unsigned) distance, Lower and Upper are
     * the limits, Weight is the weight. Returns 0.0 if Lower<=Actual<=Upper.
     * Stores the actual distance, the violated distance (Lower or Upper)
     * and the relative violation value.
     * Without arguments, returns the current relative violation value.
     */
    float rel_viol() const { return(Viol); }
    float rel_viol(float Act, float Lower, float Upper, float Weight)
    {
	Viol=0.0;
	if (Weight<=0.0 || Lower<=Actual && Actual<=Upper) return(0.0);
	Actual=Act; Strict=Weight;
	Ideal=(Actual<Lower)? Lower: Upper;
	Viol=fabsf(Ideal-Actual);
	if (Ideal>0.0) Viol/=Ideal;
	Viol*=Weight;
	return(Viol);
    }
    // END of rel_viol()

    // rel_error(): the relative error (violation divided by the weight).
    float rel_error() const { return(Viol/Strict); }
    
	// output
    friend ostream& operator<<(ostream& Out, const Viol_& V);
};
// END OF CLASS Viol_

/* Class Viollist_: stores a list of violations in descending
 * relative violation order.
 */
class Viollist_
{
    // data
    protected:
    
    List1_<Viol_> Vl;	// sorted list of violations
    
    // methods
    public:
    
	// access
    /* add_viollist(): adds the violation object V to the violation
     * list if the relative violation is larger than Minrelv. 
     * Return value: 0 if no insertion was performed, non-0 otherwise.
     */
    int add_viol(const Viol_& V, float Minrelv=0.05);
    
    /* write_file(): writes the contents of the calling object to
     * a file Outfile or to cout if Outfile==NULL or cannot be opened.
     * Return value: 1 if OK, 0 if written to cout.
     */
    int write_file(const char *Outfile) const;
    
    friend ostream& operator<<(ostream& Out, const Viollist_& Vl);
};
// END OF CLASS Viollist_

// ==== END OF HEADER Viol.h ====
#endif	/* VIOL_CLASSES */
