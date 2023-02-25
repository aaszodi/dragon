#ifndef ARRAY_TMPL_DEFS
#define ARRAY_TMPL_DEFS

// ==== FUNCTIONS Array.c++ ====

/* A general array template with optional safe indexing. */

// SGI C++ 4.0, IRIX 5.3, 7. Mar. 1996. Andris

// ==== Array_ DEFINITIONS ====

// ---- Constructors ----

/* init to N-dimensional array with undefined elements.
 * Default dim is 0.
 */
template <class T_>
Array_<T_>::Array_(unsigned int N)
{
    A=NULL; Len=0;
    if (!N) return; // empty array
    
    if ((A=new T_ [N])==NULL) prt_error(ARRAY_NO_MEM, "Array_(N)");
    else Len=N;
}

/* Init to N-dimensional array and set all elements to Item. */
template <class T_>
Array_<T_>::Array_(const T_& Item, unsigned int N)
{
    A=NULL; Len=0;
    if (!N) return; // empty array
    
    if ((A=new T_ [N])==NULL) prt_error(ARRAY_NO_MEM, "Array_(N)");
    else Len=N;
    
    for (register unsigned int i=0; i<Len; i++)
	A[i]=Item;
}

/* init to values contained in the array Arr[] which is assumed
 * to be N items long. If Arr==NULL or N==0, then a 0-dimensional
 * array is created.
 */
template <class T_>
Array_<T_>::Array_(const T_ Arr[], unsigned int N)
{
    A=NULL; Len=0; 
    if (!N || Arr==NULL) return;
    if ((A=new T_ [N])==NULL) prt_error(ARRAY_NO_MEM, "Array_(Arr, N)");
    else
    {
	Len=N; 
	if (Arr!=NULL)
	{
	    unsigned int i;
	    for (i=0; i<Len; i++) A[i]=Arr[i];
	}
    }
}

/* init by another Array_ Arr. */
template <class T_>
Array_<T_>::Array_(const Array_<T_>& Arr)
{
    A=NULL; Len=0; 
    if (!Arr.Len) return;   // init with empty array
    
    if ((A=new T_ [Arr.Len])==NULL) prt_error(ARRAY_NO_MEM, "Array_(Arr)");
    else
    {
	unsigned int i; Len=Arr.Len;
	for (i=0; i<Len; i++) A[i]=Arr.A[i];
    }
}

// ---- Safe access ----

/* Safe access is implemented by overloading the function call
 * operator (). This is a convention I'm going to stick to:
 * the index operator [] will be used for unsafe access.
 * The program is aborted on an illegal access attempt.
 */
template <class T_>
T_& Array_<T_>::operator()(unsigned int Idx)
{
    if (A==NULL) { prt_error(ARRAY_EMPTY_ACCESS, "Array(i)"); abort(); }
    if (Idx>=Len) { prt_error(ARRAY_BAD_RANGE, "Array(i)"); abort(); }
    return(A[Idx]);
}

template <class T_>
const T_& Array_<T_>::operator()(unsigned int Idx) const
{
    if (A==NULL) { prt_error(ARRAY_EMPTY_ACCESS, "Array(i)"); abort(); }
    if (Idx>=Len) { prt_error(ARRAY_BAD_RANGE, "Array(i)"); abort(); }
    return(A[Idx]);
}

// ---- Conventional (C-style) array conversion ----

/* array(N): creates an array of T_ objects and copies the
 * contents of the calling Array_ object into it. Returns the
 * address. Also returns the length of the array in N.
 * If there is no memory or the object was empty,
 * then NULL is returned and N==0.
 */
template <class T_>
T_* Array_<T_>::array(unsigned int& N) const
{
    if (!Len || A==NULL) { N=0; return(NULL); }	// calling object was empty
    
    T_ *Arr;
    unsigned int i;
    
    Arr=new T_ [Len];
    if (Arr==NULL)
    {
	prt_error(ARRAY_NO_MEM, "array(&N)");
	N=0; return(NULL);
    }
    for (i=0; i<Len; i++) Arr[i]=A[i];
    N=Len; return(Arr);
}
// END of array(N)

/* array(Arr[],N): sets the array values to those contained in Arr[]
 * which is assumed to be N items long. It is the programmer's
 * responsibility to make sure this is so. If Arr is not NULL
 * and N!=0, then the old array values are thrown away, the
 * new length will be set to N and N objects from Arr[] will be copied.
 * If either Arr==NULL or N==0 then an empty object is made.
 * Returns the calling object.
 */
template <class T_>
Array_<T_>& Array_<T_>::array(const T_ Arr[], unsigned int N)
{
    if (Arr==NULL || !N)    // make empty array
    {
	delete [] A; A=NULL; Len=0;
	return(*this);
    }
    
    if (N!=Len)
    {
	T_ *New=new T_ [N];
	if (New==NULL)	    // no mem available, do nothing
	{
	    prt_error(ARRAY_NO_MEM, "array(Arr, N)");
	    return(*this);
	}
	else
	{
	    delete [] A; A=New; Len=N;
	}
    }
    unsigned int i;
    for (i=0; i<N; i++) A[i]=Arr[i];
    return(*this);
}
// END of array(Arr[],N)

/* set_values(): sets all items to Val. Returns calling object. */
template <class T_>
Array_<T_>& Array_<T_>::set_values(const T_& Val)
{
    if (!Len || A==NULL) return(*this);	// empty array, no action is taken
    
    for (unsigned int i=0; i<Len; i++) A[i]=Val;
    return(*this);
}
// END of set_values()

// ---- Size conversion ----

/* len(): with no arguments, returns the actual length.
 * With an uint argument N, sets the length to N. If N==Len or there is
 * no memory available, no action is taken. if N<Len, then
 * only the first N items are kept. If N>Len, then the
 * items after the first Len are undefined.
 * Returns the old length.
 */
template <class T_>
unsigned int Array_<T_>::len(unsigned int N)
{
    if (N==Len) return(Len);	// same size, do nothing
    
    unsigned int Oldlen=Len;	// save old size
    if (!N)
    {	    // make empty array
	delete [] A; A=NULL; Len=0;
	return(Oldlen);
    }
    
    T_ *New;	// allocate storage for new items
    New=new T_ [N];
    if (New==NULL) { prt_error(ARRAY_NO_MEM, "len(N)"); return(Oldlen); }
    unsigned int i;
    for (i=0; i<((N>Len)? Len: N); i++) New[i]=A[i];
    delete [] A; A=New; 
    Len=N;
    return(Oldlen);
}
// END of len(N)

// ---- Assignment ----

template <class T_ >
Array_<T_>& Array_<T_>::operator=(const Array_<T_>& Arr)
{
    if (this==&Arr) return(*this);	// x=x: do nothing
    
    if (!Arr.Len)   // total destruction: x=(empty)
    {
	delete [] A; A=NULL; Len=0;
	return(*this);
    }
    
    if (Len!=Arr.Len)	// re-allocation is necessary
    {
	T_ *New=new T_ [Arr.Len];
	if (New==NULL)	    // no mem available, do nothing
	{
	    prt_error(ARRAY_NO_MEM, "=");
	    return(*this);
	}
	else
	{
	    delete [] A; A=New; Len=Arr.Len;
	}
    }
    unsigned int i;
    for (i=0; i<Len; i++) A[i]=Arr.A[i];	// copy values
    return(*this);
}
// END of operator =

// ---- Error messages ----

/* prt_error(): prints an error message to cerr. Etyp is the
 * type of the error, Funcnm is a string that contains the
 * name of the function in which the error has occurred.
 */
template <class T_>
void Array_<T_>::prt_error(Array_Errtype_ Etyp, const char *Funcnm) const
{
    cerr << "\a\n? Array_::" << Funcnm <<": ";
    switch(Etyp)
    {
	case ARRAY_NO_MEM:
	cerr << "Out of memory\n"; break;
	case ARRAY_BAD_RANGE:
	cerr << "Index out of range\n"; break;
	case ARRAY_EMPTY_ACCESS:
	cerr<< "Access to empty array attempted\n"; break;
	default:
	cerr << "Error (no code available)\n"; break;
    }
}
// END of prt_error()

// ==== END OF FUNCTIONS Array.c++ ====

#endif	/* ARRAY_TMPL_DEFS */
