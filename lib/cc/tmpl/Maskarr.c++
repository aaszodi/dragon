#ifndef MASKARR_TMPL_DEFS
#define MASKARR_TMPL_DEFS

// ==== TEMPLATE FUNCTION DEFINITIONS Maskarr.c++ ====

/* The "masked array" template class. Items are stored
 * in an array and can be switched "on" or "off". Access
 * and modification of a particular item is possible only if
 * it is switched "on" and indices are translated so that
 * [i] accesses the i-th active item.
 */

// SGI C++ 4.0, IRIX 5.3, 7. Mar. 1996. Andris Aszodi

// ==== Maskarr_ MEMBER FUNCTIONS ====

// ---- Constructors ----

/* Initialises the array to hold N items (default=0). 
 * All items will be active. This is the default ctor.
 */
template <class T_>
Maskarr_<T_>::Maskarr_(unsigned int N)
{
    if (!N)
    { Data=NULL; Idx=NULL; }
    else
    { Data=new T_ [N]; Idx=new T_* [N]; }
    Mask.len(N); Mask.set_values(true);
    update_idx();	    // assign ptrs in Idx and set Idxno
}

/* Initialises the array with Active. Maximal storage will be Active.len(),
 * the activation pattern will be Active. The array items are undefined
 * unless T_ has a default ctor.
 */
template <class T_>
Maskarr_<T_>::Maskarr_(const Bits_& Active): Mask(Active)
{
    if (Mask.len())
    { Data=new T_ [Mask.len()]; Idx=new T_* [Mask.len()]; }
    else
    { Data=NULL; Idx=NULL; }
    update_idx();
}

/* Initialises with the conventional array Arr, to have length Len.
 * It is the caller's responsibility to make sure Arr is at least Len long.
 * If Arr==NULL or Len==0 then an empty array is constructed silently.
 * All members will be active.
 */
template <class T_>
Maskarr_<T_>::Maskarr_(const T_* Arr, unsigned int Len)
{
    if (Arr==NULL || !Len)
    {
	Len=0; Data=NULL;
	Idx=NULL;
    }
    else
    {
	Data=new T_ [Len]; Idx=new T_* [Len];
	for (unsigned int i=0; i<Len; i++) Data[i]=Arr[i];  // copy
    }
    Mask.len(Len); Mask.set_values(true);
    update_idx();
}

/* The copy constructor */
template <class T_>
Maskarr_<T_>::Maskarr_(const Maskarr_<T_>& Marr)
{
    Mask=Marr.Mask;
    if (Mask.len())
    {
	Data=new T_ [Mask.len()]; Idx=new T_* [Mask.len()];
	for (unsigned int i=0; i<Mask.len(); i++)
	    Data[i]=Marr.Data[i];
    }
    else
    { Data=NULL; Idx=NULL; }
    update_idx();
}

// ---- Assignment ----

/* operator = : performs a "destructive assignment": the lvalue is
 * destroyed and an exact copy of the rvalue is built.
 */
template <class T_>
Maskarr_<T_>& Maskarr_<T_>::operator=(const Maskarr_<T_>& Marr)
{
    if (this==&Marr) return(*this);	// x=x, do nothing
    
    data_resize(Marr.Mask.len(), Marr.Data);	// copy data
    Mask=Marr.Mask; 
    update_idx();
    return(*this);
}
// END of operator =

// ---- Access ----

/* operator []: direct access to array items is not allowed.
 * [Index] will return the Index-th active item. If Index is
 * out of range or no items are active or access to empty arrays 
 * is attempted, then [] aborts with a warning.
 */
template <class T_>
const T_& Maskarr_<T_>::operator[](unsigned int Index) const
{
    if (Data==NULL || Idx==NULL)    // empty!
    {
	cerr<<"\n! Maskarr_::[]: Access to empty object attempted, abort...\n";
	abort();
    }
    if (!Idxno)
    {
	cerr<<"\n! Maskarr_["<<Index<<"]:: No active items, abort...\n";
	abort();
    }
    if (Index>=Idxno)
    {
	cerr<<"\n! Maskarr_::["<<Index<<"]: out of range [0.."<<(Idxno-1)
	    <<"], abort\n";
	abort();
    }
    return(*(Idx[Index]));
}

template <class T_>
T_& Maskarr_<T_>::operator[](unsigned int Index)
{
    if (Data==NULL || Idx==NULL)    // empty!
    {
	cerr<<"\n! Maskarr_::[]: Access to empty object attempted, abort...\n";
	abort();
    }
    if (!Idxno)
    {
	cerr<<"\n! Maskarr_["<<Index<<"]:: No active items, abort...\n";
	abort();
    }
    if (Index>=Idxno)
    {
	cerr<<"\n! Maskarr_::["<<Index<<"]: out of range [0.."<<(Idxno-1)
	    <<"], abort\n";
	abort();
    }
    return(*(Idx[Index]));
}
// END of operator []

/* active(Index): returns the activation status of the Index-th data item
 * or false if Index is out of range.
 * active(Index, Value): changes the activation status of the Index-th
 * data item. Does nothing if Index is out of range. Returns old
 * activation status.
 * NOTE: Index here indexes the Data array directly.
 */
template <class T_>
bool Maskarr_<T_>::active(unsigned int Index) const { return(Mask.get_bit(Index)); }
template <class T_>
bool Maskarr_<T_>::active(unsigned int Index, bool Value)
{
    if (!Mask.len()) return(false);	// empty: do nothing
    
    bool Oldval=Mask.set_bit(Index, Value);
    if (Index<Mask.len() && Value!=Oldval)  // valid Index, status changed
	update_idx();	// update index ptrs
    return(Oldval);
}
//END of active()

/* mask(): returns the current activation mask.
 * mask(Newmask): changes the activation mask to Newmask. If Newmask
 * has a different length, this will force the re-sizing of the Data array as well.
 * Returns old activation mask.
 * mask(Value): sets all activation bits to Value. Returns old mask.
 */
template <class T_>
Bits_ Maskarr_<T_>::mask(const Bits_& Newmask)
{
    if (Mask==Newmask) return(Mask);	// phew, don't have to do anything
    
    Bits_ Oldmask=Mask; 
    if (Newmask.len()!=Oldmask.len())
	data_resize(Newmask.len());
    Mask=Newmask; update_idx();
    return(Oldmask);
}
template <class T_>
Bits_ Maskarr_<T_>::mask(bool Value)
{
    Bits_ Oldmask=Mask; Mask.set_values(Value);
    if (Oldmask!=Mask) update_idx();
    return(Oldmask);
}
// END of mask()

// ---- Size ----

/* len(): returns the length of the Data vector.
 * len(Newlen): adjusts the size of the Data vector and the Mask. Does
 * nothing if the new length is the same as the old. If the length grows, 
 * then the new items will have undefined values and will be switched off.
 * Returns old length.
 */
template <class T_>
unsigned int Maskarr_<T_>::len(unsigned int Newlen)
{
    unsigned int Oldlen=Mask.len();
    if (Newlen==Oldlen) return(Oldlen);  // do nothing
    
    data_resize(Newlen);    // change Data size
    Mask.len(Newlen);	// change Mask size
    update_idx();  
    return(Oldlen);
}
// END of len()

// ---- Memory management ----

/* update_idx(): updates the index pointer array for the calling
 * object. The ii-th element of Idx will point to the ii-th active item in Data.
 * Must be called every time Mask or Data is modified. 
 * Return value: the no. of active items (Idxno). Protected
 */
template <class T_>
unsigned int Maskarr_<T_>::update_idx()
{
    if (Idx==NULL) return(Idxno=0);	// empty object
    
    Idxno=Mask.on_no();	// save no. of active items
    if (!Idxno) return(0);    // no active items
    
    unsigned int i, ii;
    for (i=ii=0; i<Mask.len(); i++)
	if (Mask.get_bit(i))	// active
	    Idx[ii++]=Data+i;	// store a ptr to the ii-th active item
    return(Idxno);
}
// END of update_idx()

/* data_resize(): resizes the Data and Idx arrays of the calling object to Len.
 * Does nothing if Len matches the previous length (since this is
 * obtained from Mask.len(), data_resize() MUST be called prior
 * to any modifications to Mask!). If Newarr!=NULL, then it is 
 * assumed to be a Len-long array and its contents are copied into
 * the resized Data array. If Newarr==NULL (the default), then
 * the old values are retained in Data (truncation loses some from the
 * end, lengthening appends undefined values). update_idx() must always
 * be called after this method to ensure consistency.
 * Return value: the new length. Protected
 */
template <class T_>
unsigned int Maskarr_<T_>::data_resize(unsigned int Len, const T_* Newarr)
{
    if (!Len)	// make it empty
    {
	delete [] Data; delete [] Idx;
	return(0);
    }
    
    unsigned int i, Oldlen=Mask.len();
    if (Len==Oldlen)	// no change in length
    {
	if (Newarr!=NULL)
	    for (i=0; i<Len; i++) Data[i]=Newarr[i];
	return(Len);
    }
    
    T_ *Newdat=new T_ [Len];	// allocate new size
    if (Newarr!=NULL)
	for (i=0; i<Len; i++) Newdat[i]=Newarr[i];
    else
	for (i=0; i<((Len<Oldlen)? Len: Oldlen); i++) Newdat[i]=Data[i];
    delete [] Data;
    Data=Newdat;
    
    if (Idx!=NULL) delete [] Idx;   // remove old index array
    Idx=new T_* [Len];    // allocate new index array
    
    return(Len);
}
// END of data_resize()

// ==== END OF TEMPLATE FUNCTIONS Maskarr.c++ ====

#endif	/* MASKARR_TMPL_DEFS */
