#ifndef SMARTPTR_TMPL_DECLS
#define SMARTPTR_TMPL_DECLS

// ==== TEMPLATE HEADER Smartptr.h ====

/* Implements the Smartptr_ container class which can store and
 * retrieve all derived classes of a given base class (with some
 * restrictions).
 */

// SGI C++ 4.0, IRIX 5.3, 12. Sep. 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ==== CLASSES ====

/* Class Smartptr_ : stores a dynamically allocated object derived
 * from the base class B_. Access to the object inside Smartptr_<B_>
 * is through a base class pointer.
 * IMPORTANT NOTE: Smartptr_ cannot store just anything. The trick
 * relies on "virtual constructors" (v_ctor), which are not supported in C++
 * but can be emulated. The base class B_ must declare the following
 * public pure virtual function by overloading the function call
 * operator:-
 * 
 * virtual void B_::operator()(B_*& Bptr) const =0;
 *
 * In turn, each derived class D_ should define a corresponding function:-
 * 
 *  void D_::operator()(B_*& Bptr) const
 *  {
 *	return(Bptr=new D_(*this));	// allocate an exact copy of the calling object
 *  }				// assumes that there's a D_(const D_&) ctor
 *
 * This "virtual constructor" creates an exact copy of the calling object
 * with the proper memory layout, then returns a ptr to it, converted
 * to a base class ptr. Smartptr_ keeps track of the object inside it
 * through this ptr. Smartptr_ objects must be initialised either with
 * a const D_ object or by another Smartptr_ object. A default ctor is
 * also provided but it inits the internal ptr to NULL. (Dangerous.)
 * The funny signature of the v_ctor is necessary for
 * classes with more ancestors: if D_ is derived from B1_ and B2_
 * and both base classes declare v_ctors this way, then D_ objects can be
 * stored in Smartptr_<B1_> or Smartptr_<B2_>. Of course
 * D_ must have two versions of the v_ctor:-
 *
 *  void D_::operator()(B1_*& B1ptr) const;
 *  void D_::operator()(B2_*& B2ptr) const;
 *
 * and the appropriate base class is deduced from the argument type.
 * When a Smartptr_ object dies, it destroys the member
 * object: it is vital, therefore, that ~B_() is virtual.
 * The whole point is that an Array_< Smartptr_<B_> > can store a
 * "mixed array" of D1_,  D2_, ..., Dn_ objects if they were all
 * derived from B_ and Arr[i]()->fn(...) always calls the proper virt fn.
 */
template <class B_>
class Smartptr_
{
    // data
    B_ *Ptr;	// base-class ptr to family of derived-class objects
    
    // methods
    public:
    
	// constructors
    /* Empty default constructor; dangerous but necessary since
     * arrays use this
     */
    Smartptr_(): Ptr(NULL) {}
    
    /* Init by a base class object or any of the derived class objects */
    Smartptr_(const B_& B) { B(Ptr); }
    
    /* Copy constructor */
    Smartptr_(const Smartptr_<B_>& S)
    { if (S.Ptr==NULL) Ptr=NULL; else (*S.Ptr)(Ptr); }
    
	// destructor
    ~Smartptr_() { delete Ptr; }
    
	// assignment
    Smartptr_<B_>& operator=(const Smartptr_<B_>& S);
    
	// access
    /* The pointer-to-member operator -> and the address-of operator
     * & are overloaded for accessing the
     * object within Smartptr_ through a base class pointer B_*.
     * Caution is recommended.
     */
    const B_ *operator->() const { return(Ptr); }
    B_ *operator->() { return(Ptr); }
    const B_ *operator&() const { return(Ptr); }
    B_ *operator&() { return(Ptr); }
    
    /* Dereferencing: if the object P was empty (P.Ptr==NULL), then
     * the program is aborted with an error message.
     */
    const B_& operator*() const;
    B_& operator*();
};
// END OF CLASS Smartptr_

#ifdef INCLUDE_TMPL_DEFS
#include "Smartptr.c++"
#endif

// ==== END OF TEMPLATE HEADER Smartptr.h ====

#endif	/* SMARTPTR_TMPL_DECLS */
