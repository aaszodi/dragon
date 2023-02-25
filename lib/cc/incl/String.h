#ifndef STRING_CLASS
#define STRING_CLASS

// ==== HEADER String.h ====

/* Implements a String_ class which should be the part of
 * the standard C++ library... only it isn't
 */

// SGI C++, IRIX 6.2, 29. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <ctype.h>

// ==== CLASSES ==== 

/* String_: a proper string class that has the functionality of the
 * standard ANSI C string library in <string.h>. Extra syntactic
 * sugar thrown in assignment, concatenation and equality.
 * Note that some of the fancy string functions are missing.
 */
class String_
{
    // data
    protected:
    char *S;	// \0-terminated classic C string
    unsigned int Len;	// the allocated length of S (>=1).
    
    // methods
    public:
    
	// constructors
    /* Inits to store N characters (default 0). Allocates space for N+1
     * chars, terminates with \0. Clears to all-zero.
     */
    String_(unsigned int N=0);
    
    /* Inits with a conventional C string Sc. If Sc==NULL then a 1-char
     * empty string is allocated within.
     */
    String_(const char *Sc);
    
    /* The copy constructor */
    String_(const String_& Str);
    
	// destructor
    ~String_() { delete [] S; }
    
	// assignment
    String_& operator=(const String_& Str);
    
	// conversions
    operator const char* () const { return(S); }
    
    /* long_str(): converts the (signed) long value L to its base-10
     * string representation and stores it in the calling object.
     */
    void long_str(long L);
    
	// access
    /* []: accesses the Idx-th character in the string. If Idx is
     * out of range then a warning is printed to stderr and the 0-th
     * char is returned. There is no unsafe access version.
     * Note that it is possible to access chars beyond the logical end
     * of the string (the first \0) but within the allocated length.
     */
    char operator[](int Idx) const;
    char& operator[](int Idx);
    
    /* strchr(),strrchr(): returns the index of the first/last occurrence
     * of char C in the calling object or -1 if not found. Note that this
     * behaviour is different from the C functions which return pointers.
     */
    int strchr(int C) const;
    int strrchr(int C) const;
    
    /* strstr(): returns the index of the first occurrence of Str in the
     * calling object or -1 if it was not there. Note the difference between
     * this and the ANSI C behaviour which returns a pointer.
     */
    int strstr(const String_& Str) const;
    
    /* max_len(): returns the max. number of chars in the string. This may be larger
     * than the logical length which can be accessed by strlen().
     * max_len(L): sets the length of the string to L and returns the old length.
     * If the new length is longer than the old, then the extra new chars
     * will be set to '\0', if shorter, then the overhang will be truncated, 
     * if L is the same as the old length then nothing happens.
     */
    unsigned int max_len() const { return(Len-1); }
    unsigned int max_len(unsigned int L);
    
    /* strlen(): returns the logical length of the string argument. */
    friend unsigned int strlen(const String_& Str) { return(::strlen(Str.S)); }
    
    /* !S: returns 1 if S is empty, 0 otherwise. */
    int operator!() const { return(S[0]=='\0'); }
    
	// character case
    /* tolower(), toupper(): convert all chars in the calling string
     * to lower- and uppercase,  respectively.
     * Return value: a ref. to the calling object.
     */
    String_& tolower()
    {
	for (register unsigned int i=0; i<strlen(S); i++)
	    S[i]=::tolower(S[i]);
	return(*this);
    }
    
    String_& toupper()
    {
	for (register unsigned int i=0; i<strlen(S); i++)
	    S[i]=::toupper(S[i]);
	return(*this);
    }
        
	// comparison
    /* == and != compare two strings and return a meaningful logical value (0|1)
     * unlike the counter-intuitive strcmp() they use.
     * <, <=, >, >= perform lexicographic comparison.
     */
    unsigned int operator==(const String_& Str) const { return((strcmp(S, Str.S))? 0: 1); }
    unsigned int operator!=(const String_& Str) const { return((strcmp(S, Str.S))? 1: 0); }
    unsigned int operator<(const String_& Str) const { return((strcmp(S, Str.S)<0)? 1: 0); }
    unsigned int operator<=(const String_& Str) const { return((strcmp(S, Str.S)<=0)? 1: 0); }
    unsigned int operator>(const String_& Str) const { return((strcmp(S, Str.S)>0)? 1: 0); }
    unsigned int operator>=(const String_& Str) const { return((strcmp(S, Str.S)>=0)? 1: 0); }
    
	// concatenation
    /* S+=Str: appends Str to S and returns S.
     * S+Str: appends Str to S and returns the resulting object.
     */
    String_& operator+=(const String_& Str);
    String_ operator+(const String_& Str) const;

	// input/output
    /* In>>Str: Skips whitespaces and reads non-whitespace chars into
     * Str until it gets filled up or a whitespace is encountered. Will NOT
     * adjust the length of Str dynamically! 
     */
    friend istream& operator>>(istream& In, String_& Str);
    
    /* Out<<Str: writes the string Str to Out. */
    friend ostream& operator<<(ostream& Out, const String_& Str);
};
// END OF CLASS String_

// ==== END OF HEADER String.h ====

#endif	/* STRING_CLASS */
