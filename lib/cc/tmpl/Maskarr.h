#ifndef MASKARR_TMPL_DECLS
#define MASKARR_TMPL_DECLS

// ==== TEMPLATE HEADER Maskarr.h ====

/* The "masked array" template class. Items are stored
 * in an array and can be switched "on" or "off". Access
 * and modification of a particular item is possible only if
 * it is switched "on" and indices are translated so that
 * [i] accesses the i-th active item.
 */

// SGI C++ 4.0, IRIX 5.3, 7. Mar. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- INCLUDE FILES ----

#include "Bits.h"   // the bit vector class

// ==== CLASSES ==== 

/* Class Maskarr_ : the masked array template class.
 * The elements of the array can be active or inactive;
 * activation status is held in a bit vector. Only active
 * items can be accessed or modified.
 */
template <class T_>
class Maskarr_
{
    // data
    protected:
    
    T_ *Data;	// the array of the items (length is Mask.len())
    T_ **Idx;	// array of pointers to active items, NULL if all inactive
    unsigned int Idxno;	    // no. of active items ( same as Mask.on_no() )
    Bits_ Mask;	// the bit-vector holding the activation status
    
    // methods
    public:
    
	// constructors
    /* Initialises the array to hold N items (default=0). 
     * All items will be active. This is the default ctor.
     */
    Maskarr_(unsigned int N=0);
    
    /* Initialises the array with Active. Maximal storage will be Active.len(),
     * the activation pattern will be Active. The array items are undefined
     * unless T_ has a default ctor.
     */
    Maskarr_(const Bits_& Active);
    
    /* Initialises with the conventional array Arr, to have length Len.
     * It is the caller's responsibility to make sure Arr is at least Len long.
     * If Arr==NULL or Len==0 then an empty array is constructed silently.
     * All members will be active.
     */
    Maskarr_(const T_* Arr, unsigned int Len);
    
    /* The copy constructor */
    Maskarr_(const Maskarr_<T_>& Marr);
    
	// destructor
    virtual ~Maskarr_() { delete [] Data; delete [] Idx; }
    
	// assignment 
    /* operator = : performs a "destructive assignment": the lvalue is
     * destroyed and an exact copy of the rvalue is built.
     */
    Maskarr_<T_>& operator=(const Maskarr_<T_>& Marr);
    
	// access
    /* operator []: direct access to array items is not allowed.
     * [Index] will return the Index-th active item. If Index is
     * out of range or no items are active or access to empty arrays 
     * is attempted, then [] aborts with a warning.
     */
    const T_& operator[](unsigned int Index) const;
    T_& operator[](unsigned int Index);
	
    /* active(Index): returns the activation status of the Index-th data item
     * or false if Index is out of range.
     * active(Index, Value): changes the activation status of the Index-th
     * data item. Does nothing if Index is out of range. Returns old
     * activation status.
     * NOTE: Index here indexes the Data array directly.
     */
    bool active(unsigned int Index) const;
    bool active(unsigned int Index, bool Value);
    
    /* mask(): returns the current activation mask.
     * mask(Newmask): changes the activation mask to Newmask. If Newmask
     * has a different length, this will force the re-sizing of the Data array as well.
     * Returns old activation mask.
     * mask(Value): sets all activation bits to Value. Returns old mask.
     */
    const Bits_& mask() const { return(Mask); }
    Bits_ mask(const Bits_& Newmask);
    Bits_ mask(bool Value);
    
	// size
    
    /* len(): returns the length of the Data vector.
     * len(Newlen): adjusts the size of the Data vector and the Mask. Does
     * nothing if the new length is the same as the old. If the length grows, 
     * then the new items will have undefined values and will be switched off.
     * Returns old length.
     */
    unsigned int len() const { return(Mask.len()); }
    unsigned int len(unsigned int Newlen);
    
    /* active_len(): returns the number of active elements. */
    unsigned int active_len() const { return(Idxno); }

    // protected methods
    protected:
    
    unsigned int update_idx();
    unsigned int data_resize(unsigned int Len, const T_* Newarr=NULL);

};
// END OF CLASS Maskarr_

#ifdef INCLUDE_TMPL_DEFS
#include "Maskarr.c++"
#endif

// ==== END OF TEMPLATE HEADER Maskarr.h ====

#endif	/* MASKARR_TMPL_DECLS */
