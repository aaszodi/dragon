#ifndef SMARTPTR_TMPL_DEFS
#define SMARTPTR_TMPL_DEFS

// ==== TEMPLATE METHODS Smartptr.c++ ====

/* Implements the Smartptr_ container class which can store and
 * retrieve all derived classes of a given base class (with some
 * restrictions).
 */

// SGI C++ 4.0, IRIX 5.3, 27. Apr. 1995. Andris Aszodi

// ==== Smartptr_ METHODS ====

// ---- Assignment ----

template <class B_>
Smartptr_<B_>& Smartptr_<B_>::operator=(const Smartptr_<B_>& S)
{
    if (Ptr==S.Ptr) return(*this);    // x=x
    delete Ptr; 
    if (S.Ptr==NULL) Ptr=NULL;	// copy empty object
    else (*S.Ptr)(Ptr);	// v_ctor call
    return(*this);
}
// END of =

// ---- Access ----

/* Dereferencing: if the object P was empty (P.Ptr==NULL), then
 * the program is aborted with an error message. 
 */
template <class B_>
const B_& Smartptr_<B_>::operator*() const
{
    if (Ptr==NULL)  // problem
    {
	cerr<<"! *Smartptr_<_>:: Accessing empty const object\n";
	abort();
    }
    return(*Ptr);
}

template <class B_>
B_& Smartptr_<B_>::operator*()
{
    if (Ptr==NULL)  // big problem
    {
	cerr<<"! *Smartptr_<_>:: Accessing empty object, prepare for disaster!\n";
	abort();
    }
    return(*Ptr);
}
// END of *()

// ==== END OF METHODS Smartptr.c++ ====

#endif	/* SMARTPTR_TMPL_DEFS */
