#ifndef ARRAY_TMPL_DECLS
#define ARRAY_TMPL_DECLS

// ==== HEADER Array.h ====

/* A general array template with optional safe indexing. */

// SGI C++ 4.0, IRIX 5.3, 7. Mar. 1996. Andris

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>

// ==== CLASSES ====

/* Class Array_: a simple 1D array template with safe
 * indexing and re-allocation capabilities.
 */
template <class T_>
class Array_
{
    // data
    T_ *A;	// address of the array (NULL if empty)
    unsigned int Len; // length of array (0 if empty)
    
    // functions
    public:
    
	// constructors
    /* init to N-dimensional array with undefined elements.
     * Default dim is 0.
     */
    Array_(unsigned int N=0);
    
    /* Init to N-dimensional array and set all elements to Item. */
    Array_(const T_& Item, unsigned int N);
    
    /* init to values contained in the array Arr[] which is assumed
     * to be N items long. If Arr==NULL or N==0, then a 0-dimensional
     * array is created.
     */
    Array_(const T_ Arr[], unsigned int N);
    
    /* init by another Array_ Arr. */
    Array_(const Array_<T_>& Arr);
    
	// destructor
    ~Array_() { delete [] A; Len=0; }
    
	// unsafe member access
    /* The Idx value is not checked. */
    T_& operator[](unsigned int Idx) { return(A[Idx]); }
    const T_& operator[](unsigned int Idx) const { return(A[Idx]); };
    
	// safe member access
    
    /* Safe access is implemented by overloading the function call
     * operator (). This is a convention I'm going to stick to:
     * the index operator [] will be used for unsafe access.
     */
    T_& operator()(unsigned int Idx);
    const T_& operator()(unsigned int Idx) const;
    
	// conventional (C-style) array conversion
    /* array(N): creates an array of T_ objects and copies the
     * contents of the calling Array_ object into it. Returns the
     * address. Also returns the length of the array in N.
     * If there is no memory or the object was empty,
     * then NULL is returned and N==0.
     */
    T_* array(unsigned int& N) const;
    
    /* array(Arr[],N): sets the array values to those contained in Arr[]
     * which is assumed to be N items long. It is the programmer's
     * responsibility to make sure this is so. If Arr is not NULL
     * and N!=0, then the old array values are thrown away, the
     * new length will be set to N and N objects from Arr[] will be copied.
     * If either Arr==NULL or N==0 then an empty object is made.
     * Returns the calling object.
     */
    Array_<T_>& array(const T_ Arr[], unsigned int N);
    
    /* set_values(): sets all items to Val. Returns calling object. */
    Array_<T_>& set_values(const T_& Val);
    
	// size access
    /* len(): with no arguments, returns the actual length.
     * With an uint argument N, sets the length to N. If N==Len or there is
     * no memory available, no action is taken. if N<Len, then
     * only the first N items are kept. If N>Len, then the
     * items after the first Len are undefined.
     * Returns the old length.
     */
    unsigned int len() const { return(Len); }
    unsigned int len(unsigned int N);

	// assignment
    Array_<T_>& operator=(const Array_<T_>& Arr);

    private:
	// error messages
    enum Array_Errtype_ {ARRAY_NO_MEM, ARRAY_BAD_RANGE, ARRAY_EMPTY_ACCESS};
    void prt_error(Array_Errtype_ Etyp, const char *Funcnm) const;
    
};
// END OF CLASS Array_

#ifdef INCLUDE_TMPL_DEFS
#include "Array.c++"
#endif

// ==== END OF HEADER Array.h ====

#endif	/* ARRAY_TMPL_DECLS */
