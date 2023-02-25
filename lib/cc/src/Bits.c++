// ==== FUNCTIONS Bits.c++ ====

/* A simple bit-array class for storing Boolean data. */

// SGI C++ 4.0, IRIX 5.3, 28. Sept. 1995. Andris 

// ---- HEADER ----

#include "Bits.h"

// ==== Bits_ MEMBER FUNCTIONS ====

// ---- Memory handling ----

/* alloc_arr(): creates an uint array containing Size items
 * and returns a ptr to it (or NULL if allocation failed).
 * If Arr!=NULL, then its contents (up to Size uints) are
 * copied into the new array. Arr==NULL by default. Private static
 */
inline
unsigned int* Bits_::alloc_arr(unsigned int Size, const unsigned int *Arr)
{
    if (!Size) return(NULL);	// do nothing
    
    unsigned int *New=new unsigned int [Size];
    if (New==NULL)
    {	// should throw exception here
	cerr<<"\n! Bits_::alloc_arr(): Out of memory\n";
	return(NULL);
    }
    if (Arr!=NULL) memcpy(New, Arr, Size*sizeof(unsigned int));
    return(New);
}
// END of alloc_arr()

/* zero_mask(): the last uint in the B array may be incompletely filled.
 * This method returns an uint which has 1-s for the valid positions
 * and 0-s for the unused tail. For zeroing that tail, just & to the
 * last uint. Private
 */
inline
unsigned int Bits_::zero_mask() const
{
    register unsigned int Lastbits=Bs % UINT_SIZE,  // no. of bits in last active uint
	Zmask=~0;  // all 111...1
    if (Lastbits)
    {
	// make 00..011..11 when not full
	if (Lastbits<=UINT_SIZE/2) { Zmask<<=Lastbits; Zmask=~Zmask; }
	else Zmask>>=UINT_SIZE-Lastbits;
    }
    return(Zmask);  // all 1 if Lastbits==0
}
// END of zero_mask()

// ---- Constructors ----

/* Init to hold N bits (default 0). Set all bits to 
 * Val (false by default).
 */
Bits_::Bits_(unsigned int N, bool Val):
    Bs(N), Cs(uints_needed(N))
{
    if (!N) { B=NULL; Cs=0; return; }    // empty
    
    B=alloc_arr(Cs);
    
    // this sets the extra bits at the upper end of the last uint, too!
    if (B!=NULL)
	memset(B, (!Val)? '\0': ~'\0', Cs*sizeof(unsigned int));
}

/* Init by another Bits_ object Bits (copy constructor). */
Bits_::Bits_(const Bits_& Bits):
    B(alloc_arr(Bits.Cs, Bits.B)), Bs(Bits.Bs), Cs(Bits.Cs)
{
    if (!Bits.Bs || B==NULL) return;    // empty
}

/* Init by a string Str which should have the format "110101001",
 * 0 for "false", 1 for "true". Any other char elicits a warning
 * and will be interpreted as "true". The null string "" is accepted
 * as an empty bitvector.
 */
Bits_::Bits_(const char *Str):
    Bs(strlen(Str))
{
    if (!Bs) { B=NULL; Cs=0; return; }	// empty
    
    Cs=uints_needed(Bs);
    B=alloc_arr(Cs);
    if (B!=NULL)
    {
	memset(B, 0, Cs*sizeof(unsigned int));	// zero all
	for (register int i=0; i<Bs; i++)
	{
	    char Ch=Str[Bs-i-1];
	    if (Ch=='0') continue;
	    if (Ch!='1')
		cerr<<"\n? Bits_(str): Invalid char \'"<<Ch<<"\', \'1\' used\n";
	    set_bit(i);
	}
    }
}

// ---- Assignment ----

/* Performs 'destructive assignment', i.e. the original object
 * is destroyed and then rebuilt as an exact copy of the Other object.
 */
Bits_& Bits_::operator=(const Bits_& Other)
{
    if (this==&Other) return(*this);	// x==x, so x=x is pointless
    
    if (Cs<Other.Cs)	// realloc is necessary
    {
	unsigned int *Bnew=alloc_arr(Other.Cs, Other.B);
	if (Bnew==NULL) return(*this);
	delete [] B; B=Bnew; Cs=Other.Cs;
    }
    else memcpy(B, Other.B, Other.Cs*sizeof(unsigned int));
    
    Bs=Other.Bs;
    return(*this);
}
// END of =

// ---- Equality ----

/* ==,!=: return an appropriate Boolean value. Two Bits_ arrays are
 * equal if they have the same number of bits and all bits are set
 * in the same way.
 */
bool Bits_::operator==(const Bits_& Other) const
{
    if (this==&Other) return(true);	// x==x
    if (Bs!=Other.Bs) return(false);  // different sizes
    
    // check the "full" uint's in the array first
    register unsigned int i, Uino=uints_needed(Bs);
    for (i=0; i<Uino-1; i++)
	if (B[i]!=Other.B[i]) return(false);
    
    // check the bits in the last uint
    for (i=(Uino-1)*sizeof(unsigned int); i<Bs; i++)
	if (get_bit(i)!=Other.get_bit(i)) return(false);
    return(true);
}
// END of ==

bool Bits_::operator!=(const Bits_& Other) const
{
    if (this==&Other) return(false);	// x!=x
    if (Bs!=Other.Bs) return(true);  // different sizes
    
    // check the "full" uint's in the array first
    register unsigned int i, Uino=uints_needed(Bs);
    for (i=0; i<Uino-1; i++)
	if (B[i]!=Other.B[i]) return(true);
    
    // check the bits in the last uint
    for (i=(Uino-1)*sizeof(unsigned int); i<Bs; i++)
	if (get_bit(i)!=Other.get_bit(i)) return(true);
    return(false);
}
// END of !=

// ---- Access and size ----

/* get_mask(): auxiliary function for "locating" individual bits.
 * The Idx-th bit is accessed so that the offset of the unsigned int
 * array is returned in Offs and an unsigned int bitmask is generated
 * which has 1 at the corresponding position and 0 bits otherwise.
 * NOTE: Idx is not checked here (assumed to be correct!).
 * Return value: the bitmask or 0 if Idx was out of range. Private static
 */
inline
unsigned int Bits_::get_mask(unsigned int Idx, unsigned int& Offs)
{
    Offs=Idx/UINT_SIZE;
    register unsigned int Mask=1; Mask<<=(Idx % UINT_SIZE);
    return(Mask);
}
// END of get_mask()

/* get_bit(): returns the Boolean value of the Idx-th bit.
 * If Idx is out of range, false will be returned.
 */
bool Bits_::get_bit(unsigned int Idx) const
{
    // return 0 silently if Idx is out of range or calling object is empty
    if (!Bs || Idx>=Bs) return(false);
    
    register unsigned int Offs, Mask=get_mask(Idx, Offs);
    return(bool((Mask & B[Offs])!=0));    // to convert to 0|1
}
// END of get_bit()

/* set_bit(): sets the Idx-th bit to Value (default is true).
 * If Idx is out of range, no action
 * will be taken. Returns old value or false if out-of-range.
 */
bool Bits_::set_bit(unsigned int Idx, bool Value)
{
    // return 0 silently if Idx is out of range or calling object is empty
    if (!Bs || Idx>=Bs) return(false);
    
    register unsigned int Offs, Mask=get_mask(Idx, Offs);
    bool Oldbit=bool((Mask & B[Offs])!=0);
    if (Value) B[Offs]|=Mask; else B[Offs]&=~Mask;
    return(Oldbit);
}
// END of set_bit()

/* len(Len): adjusts the size of the bit array to Len.
 * If Len==Bs,  no action is taken, if Len<Bs then the
 * tail of the array is lost, if Len>Bs then the tail of the
 * new array will be filled with false values (0).
 * Returns the old length.
 */
unsigned int Bits_::len(unsigned int Len)
{
    if (Len==Bs) return(Bs);    // no change
    
    unsigned int Oldlen=Bs;
    
    if (Len>Bs)	    // bit size grows (may fit into original Cs-long array)
    {
	register unsigned int Csnew=uints_needed(Len),  // new uint size
	    Csmin=uints_needed(Bs);	// minimal uint size for old
	
	if (Cs<Csnew)  // realloc is necessary (uint size grows as well)
	{
	    unsigned int *Bnew=alloc_arr(Csnew);	// create new array
	    if (Bnew==NULL) return(Oldlen);
	    
	    if (B!=NULL)	// there was a previous array (i.e. Bs>0,Cs>0)
	    {
		memcpy(Bnew, B, Csmin*sizeof(unsigned int));	// copy old
		delete [] B;
	    }
	    B=Bnew;	// replace old array w/ new
	}
	memset(B+Csmin, 0, (Csnew-Csmin)*sizeof(unsigned int));   // zero new uints
	
	// zero the last bits in the last (old) uint if necessary
	if (Csmin) B[Csmin-1]&=zero_mask();
	
	if (Cs<Csnew) Cs=Csnew;	// Cs had to be kept until now, cf zeroing above
    }
    
    Bs=Len;
    return(Oldlen);
}
// END of len(Len)

/* set_values(): sets all bits to Value (default is true). */
void Bits_::set_values(bool Value)
{
    if (!Bs) return;
    unsigned char Vals=(Value)? ~0: 0;
    memset(B, Vals, uints_needed(Bs)*sizeof(unsigned int));
}

/* on_no(), off_no(): return the number of bits that are ON or OFF,
 * respectively.
 */
unsigned int Bits_::on_no() const
{
    unsigned int No, i;
    for (No=0, i=0; i<Bs; i++)
	if (get_bit(i)) No++;
    return(No);
}
// END of on_no()

unsigned int Bits_::off_no() const
{
    unsigned int No, i;
    for (No=0, i=0; i<Bs; i++)
	if (!get_bit(i)) No++;
    return(No);
}
// END of off_no()

// ---- Bitwise operations ----

/* The following overlaid operators (~,&,|,^) perform exactly the same
 * bitwise operations as their standard counterparts. The bit array sizes
 * must be the same otherwise a dim mismatch error occurs, a
 * warning is printed and the left operand is returned.
 */

Bits_& Bits_::operator~()
{
    for (unsigned int i=0; i<uints_needed(Bs); i++) B[i]=~B[i];
    return(*this);
}

Bits_& Bits_::operator&=(const Bits_& Bits)
{
    if (Bs!=Bits.Bs) cerr<<"\n? Bits_::x&=y: dim mismatch\n";
    else for (unsigned int i=0; i<uints_needed(Bs); i++) B[i]&=Bits.B[i];
    return(*this);
}

Bits_& Bits_::operator|=(const Bits_& Bits)
{
    if (Bs!=Bits.Bs) cerr<<"\n? Bits_::x|=y: dim mismatch\n";
    else for (unsigned int i=0; i<uints_needed(Bs); i++) B[i]|=Bits.B[i];
    return(*this);
}

Bits_& Bits_::operator^=(const Bits_& Bits)
{
    if (Bs!=Bits.Bs) cerr<<"\n? Bits_::x^=y: dim mismatch\n";
    else for (unsigned int i=0; i<uints_needed(Bs); i++) B[i]^=Bits.B[i];
    return(*this);
}

Bits_ Bits_::operator&(const Bits_& Bits) const
{
    if (Bs!=Bits.Bs) { cerr<<"\n? Bits_::x&y: dim mismatch\n"; return(*this); }
    Bits_ Temp=(*this); Temp&=Bits; return(Temp);
}

Bits_ Bits_::operator|(const Bits_& Bits) const
{
    if (Bs!=Bits.Bs) { cerr<<"\n? Bits_::x|y: dim mismatch\n"; return(*this); }
    Bits_ Temp=(*this); Temp|=Bits; return(Temp);
}

Bits_ Bits_::operator^(const Bits_& Bits) const
{
    if (Bs!=Bits.Bs) { cerr<<"\n? Bits_::x^y: dim mismatch\n"; return(*this); }
    Bits_ Temp=(*this); Temp^=Bits; return(Temp);
}

// ---- Shifts ----

/* The following operators (<<,<<=,>>,>>=) do exactly the same things as
 * their original C counterparts. However, the rhs-s (although the type
 * is int) are not allowed to have negative values: no action is taken and
 * a warning is printed in these cases. For << and <<=, 0 bits are shifted 
 * in at the right; for >> and >>=, 0 bits are shifted in from the left.
 * Empty bit-vectors are left unchanged.
 */

Bits_& Bits_::operator<<=(int Shift)
{
    if (!Bs || !Shift) return(*this);	// empty or no shift
    if (Shift<0)
    {
	cerr<<"\n? Bits_<<="<<Shift<<": Negative shift\n";
	return(*this);
    }

    // if Shift>=UINT_SIZE, full uint's must be shifted first
    register unsigned int i, Csmin=uints_needed(Bs), Ush=Shift/UINT_SIZE;
    
    if (Ush)	// do uint-left-shift
    {
	for (i=Csmin-1; i>=Ush; i--) B[i]=B[i-Ush];
	memset(B, 0, Ush*sizeof(unsigned int));	    // zero the Ush rightmost uints
    }

    // do bitwise shifts, carrying over bits at each left side
    register unsigned int Bsh=Shift % UINT_SIZE;
    if (Bsh)
    {
	register unsigned int Carry, Prevcarry=0, Cmask=~0;	// all 11..1
	Cmask>>=Bsh; Cmask=~Cmask;  // 11..10000..0
	for (i=Ush; i<Csmin-Ush; i++)
	{
	    Carry=B[i] & Cmask;	// save leftmost Bsh bits
	    B[i] <<= Bsh;   // do the shift, rightmost Bsh bits are 0
	    B[i] |= Prevcarry;	// put there saved bits from previous shift
	    Prevcarry=Carry;	// save for next
	    Prevcarry>>=(UINT_SIZE-Bsh);    // move from left to right
	}
    }
    return(*this);
}
// END of <<=

Bits_& Bits_::operator>>=(int Shift)
{
    if (!Bs || !Shift) return(*this);	// empty or no shift
    if (Shift<0)
    {
	cerr<<"\n? Bits_>>="<<Shift<<": Negative shift\n";
	return(*this);
    }

    // if Shift>=UINT_SIZE, full uint's must be shifted first
    register unsigned int Csmin=uints_needed(Bs), Ush;
    register int i;
    
    /* subtract the "tail" Csmin % Bs bits from Shift and do
     * uint-right-shifts only if Shift is larger than this
     */
    Ush=(Csmin % Bs >= Shift)? 0: (Shift-Csmin % Bs)/UINT_SIZE;
    if (Ush)	// do uint-right-shift
    {
	for (i=0; i<Csmin-Ush; i++) B[i]=B[i+Ush];
	memset(B+Csmin-Ush, 0, Ush*sizeof(unsigned int));	    // zero the Ush leftmost uints
    }

    // do bitwise shifts, carrying over bits at each right side
    register unsigned int Bsh=Shift % UINT_SIZE;
    if (Bsh)
    {
	register unsigned int Carry, Prevcarry=0, Cmask=~0;	// all 11..1
	Cmask<<=Bsh; Cmask=~Cmask;  // 00..0011..1111
	B[Csmin-1]&=zero_mask();    // zero tail
	for (i=Csmin-Ush-1; i>=0; i--)    // proper -1 check needed, i not uint
	{
	    Carry=B[i] & Cmask;	// save rightmost Bsh bits
	    B[i] >>= Bsh;   // do the shift, leftmost Bsh bits are 0
	    B[i] |= Prevcarry;	// put there saved bits from previous shift
	    Prevcarry=Carry;	// save for next
	    Prevcarry<<=(UINT_SIZE-Bsh);    // move from right to left
	}
    }
    return(*this);
}
// END of >>=

Bits_ Bits_::operator<<(int Shift) const
{
    Bits_ Temp(*this); Temp<<=Shift;
    return(Temp);
}
// END of <<

Bits_ Bits_::operator>>(int Shift) const
{
    Bits_ Temp(*this); Temp>>=Shift;
    return(Temp);
}
// END of >>

// ---- Printing ----

/* list_bits(): lists the bits neatly to Out, UINT_SIZE bits per row.
 * Prints Fch (default '0') and Tch (default '1') for false and
 * true values,  respectively.
 */
void Bits_::list_bits(ostream& Out, char Fch, char Tch) const
{
    unsigned int i, j, Mask, Lstb=Bs % UINT_SIZE, Uino=uints_needed(Bs); 
    const unsigned int Dw=sizeof(unsigned int)*(CHAR_BIT+1);
    
    // print --- on top
    for (j=0; j<Dw; j++) Out.put('-'); Out<<endl;
    
    // print "full" uints
    for (i=0; i<((Lstb)? Uino-1: Uino); i++)
    {
	Mask=1U<<(UINT_SIZE-1);	// 1000...000
	for (j=0; j<UINT_SIZE; j++)
	{
	    Out.put((B[i] & Mask)? Tch: Fch);
	    if (!((j+1) % CHAR_BIT)) Out.put(' ');  // separate bytes
	    Mask>>=1;
	}
	Out<<'['<<(i+1)*UINT_SIZE-1<<'-'<<i*UINT_SIZE<<"]\n";
    }
    
    // print last uint (might be incomplete)
    if (Lstb)
    {
	Mask=1<<(Lstb-1);	// 1000...000 (shorter)
	for (j=0; j<UINT_SIZE; j++)
	{
	    if (j<UINT_SIZE-Lstb) Out.put(' ');	    // just left pad
	    else
	    {
		Out.put((B[Uino-1] & Mask)? Tch: Fch);
		Mask>>=1;
	    }
	    if (!((j+1) % CHAR_BIT)) Out.put(' ');
	}
	Out<<'['<<Bs-1<<'-'<<(Cs-1)*UINT_SIZE<<"]\n";
    }
    
    // print === at bottom
    for (j=0; j<Dw; j++) Out.put('='); Out<<"\n\n";
}
// END of list_bits()

/* <<: overlaid print operator that calls list_bits()
 * with the default arguments.
 */
ostream& operator<<(ostream& Out, const Bits_& Bits)
{
    Bits.list_bits(Out); return(Out);
}
// END of <<

// ==== END OF FUNCTIONS Bits.c++ ====
