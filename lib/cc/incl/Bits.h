#ifndef BITS_CLASS
#define BITS_CLASS

// ==== HEADER Bits.h ====

/* A simple bit-array class for storing Boolean data. */

// SGI C++ 7.1, IRIX 6.2, 4. Mar. 1997. Andris 

// ---- HEADERS ----

#include <stdlib.h>
#include <limits.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>

// ---- DEFINITIONS ----

#ifndef CHAR_BIT
#define CHAR_BIT 8
#endif
#define UINT_SIZE (CHAR_BIT*sizeof(unsigned int))

// ---- GLOBAL TYPES ----

/* NOTE: The SGI N32/N64 compilers recognise the built-in bool type
 * and define the macro _BOOL. This workaround is provided for
 * the old-style SGI O32 compiler.
 * Define _BOOL on the command line for 'bool'-aware non-SGI compilers.
 */
#if !(defined(_BOOL) || defined(HAS_SEEN_BOOL))
#define HAS_SEEN_BOOL
typedef unsigned char bool;
static const bool false=0;
static const bool true=1;
#endif	/* _BOOL */

// ==== CLASSES ====

/* Bits_ : A simple bit-array class. Bits are stored in 
 * an array of unsigned ints (char size CHAR_BIT is taken from
 * <limits.h>). Individual bits can be set and tested
 * and global bitwise Boolean operations can be performed.
 * Empty arrays can also be constructed.
 */
class Bits_
{
    // data
    unsigned int *B;	// the bit array
    unsigned int Cs, Bs;	// no. of uints and bits
    
    // methods
    public:
    
	// constructors
    /* Init to hold N bits (default 0). Set all bits to 
     * Val (false by default).
     */
    Bits_(unsigned int N=0, bool Val=false);
    
    /* Init by another Bits_ object Bits (copy constructor). */
    Bits_(const Bits_& Bits);
	
    /* Init by a string Str which should have the format "110101001",
     * 0 for "false", 1 for "true". Any other char elicits a warning
     * and will be interpreted as "true". The null string "" is accepted
     * as an empty bitvector.
     */
    Bits_(const char *Str);
    
	// destructor
    ~Bits_() { delete [] B; }
	
	// assignment
    /* Performs 'destructive assignment', i.e. the original object
     * is destroyed and then rebuilt as an exact copy of the Other object.
     */
    Bits_& operator=(const Bits_& Other);
    
	// equality
    /* ==,!=: return an appropriate Boolean value. Two Bits_ arrays are
     * equal if they have the same number of bits and all bits are set
     * in the same way.
     */
    bool operator==(const Bits_& Other) const;
    bool operator!=(const Bits_& Other) const;
    
	// access and size
    unsigned int cno() const { return(Cs); }
    
    /* get_bit(): returns the Boolean value of the Idx-th bit.
     * If Idx is out of range, false will be returned.
     */
    bool get_bit(unsigned int Idx) const;
    
    /* set_bit(): sets the Idx-th bit to Value (default is true).
     * If Idx is out of range, no action
     * will be taken. Returns old value or false if out of range.
     */
    bool set_bit(unsigned int Idx, bool Value=true);
    
    /* set_values(): sets all bits to Value (default is true). */
    void set_values(bool Value=true);
    
    /* len(): returns the no. of bits in the array.
     * len(Len): adjusts the size of the bit array to Len.
     * If Len==Bs,  no action is taken, if Len<Bs then the
     * tail of the array is lost, if Len>Bs then the tail of the
     * new array will be filled with false values (0).
     */
    unsigned int len() const { return(Bs); }
    unsigned int len(unsigned int Len);
    
    /* on_no(), off_no(): return the number of bits that are ON or OFF,
     * respectively.
     */
    unsigned int on_no() const;
    unsigned int off_no() const;
    
	// bitwise operations
    /* The following overlaid operators (~,&,|,^) perform exactly the same
     * bitwise operations as their standard counterparts. The bit array sizes
     * must be the same otherwise a dim mismatch error occurs, a
     * warning is printed and the left operand is returned.
     */
    Bits_& operator~();
    Bits_& operator&=(const Bits_& Bits);
    Bits_& operator|=(const Bits_& Bits);
    Bits_& operator^=(const Bits_& Bits);
    Bits_ operator&(const Bits_& Bits) const;
    Bits_ operator|(const Bits_& Bits) const;
    Bits_ operator^(const Bits_& Bits) const;

	// shifts
    /* The following operators (<<,<<=,>>,>>=) do exactly the same things as
     * their original C counterparts. However, the rhs-s (although the type
     * is int) are not allowed to have negative values: no action is taken and
     * a warning is printed in these cases. For << and <<=, 0 bits are shifted 
     * in at the right; for >> and >>=, 0 bits are shifted in from the left.
     * Empty bit-vectors are left unchanged.
     */
    Bits_& operator<<=(int Shift);
    Bits_& operator>>=(int Shift);
    Bits_ operator<<(int Shift) const;
    Bits_ operator>>(int Shift) const;

	// printing
    /* list_bits(): lists the bits neatly to Out, UINT_SIZE bits per row.
     * Prints Fch (default '0') and Tch (default '1') for false and
     * true values,  respectively.
     */
    void list_bits(ostream& Out, char Fch='0', char Tch='1') const;
    
    /* <<: overlaid print operator that calls list_bits()
     * with the default arguments.
     */
    friend ostream& operator<<(ostream& Out, const Bits_& Bits);

    // end of public methods interface
    
    private:
    static unsigned int uints_needed(unsigned int Bitno)    // minimal no. of uints
    { return(Bitno? (Bitno-1)/UINT_SIZE+1: 0); }
    static unsigned int get_mask(unsigned int Idx, unsigned int& Offs);
    
    static unsigned int* alloc_arr(unsigned int Size, const unsigned int *Arr=NULL);
    unsigned int zero_mask() const;

};
// END OF CLASS Bits_

// ==== END OF HEADER Bits.h ====

#endif  /* BITS_CLASS */
