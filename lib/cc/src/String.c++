// ==== METHODS String.c++ ====

/* Implements a String_ class which should be the part of
 * the standard C++ library... only it isn't
 */

// SGI C++, IRIX 6.2, 29. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <strstream.h>

// ---- CLASS HEADER ----

#include "String.h"

// ==== String_ METHODS ====

// ---- Constructors ----

/* Inits to store N characters (default 0). Allocates space for N+1
 * chars, terminates with \0. Clears to all-zero.
 */
String_::String_(unsigned int N): Len(N+1)
{
    S=new char [N+1]; memset(S, 0, N+1);
}

/* Inits with a conventional C string Sc. If Sc==NULL then a 1-char
 * empty string is allocated within.
 */
String_::String_(const char *Sc)
{
    if (Sc==NULL)
    {
	S=new char [Len=1]; *S='\0';
    }
    else
    {
	S=new char [Len=strlen(Sc)+1]; S[Len-1]='\0';
	strcpy(S, Sc);
    }
}

/* The copy constructor */
String_::String_(const String_& Str)
{
    S=new char [Len=Str.Len]; strcpy(S, Str.S);
    S[Len-1]='\0';
}

/* Assignment */
String_& String_::operator=(const String_& Str)
{
    if (this==&Str) return(*this);  // x=x
    if (Len!=Str.Len)
    {
	// realloc needed
	delete [] S;
	S= new char [Len=Str.Len];
    }
    memcpy(S, Str.S, Len);
    return(*this);
}
// END of =

// ---- num->str conversion ----

/* long_str(): converts the (signed) long value L to its base-10
 * string representation and stores it in the calling object.
 */
void String_::long_str(long L)
{
    ostrstream Os;	// use incore formatting
    Os<<L;  // write the value, let the system do the job
#ifdef __linux__
    Os<<'\0';
#endif
    delete [] S;    // get rid of old value
    S=Os.str();	// freeze output array and return its ptr
    Len=strlen(S)+1;	// store new maxlen
}
// END of long_str()

// ---- Access ----

/* []: accesses the Idx-th character in the string. If Idx is
 * out of range then a warning is printed to stderr and the 0-th
 * char is returned. There is no unsafe access version.
 * Note that it is possible to access chars beyond the logical end
 * of the string (the first \0) but within the allocated length.
 */
char String_::operator[](int Idx) const
{
    if (Idx<0 || Idx>=max_len())
    {
	cerr<<"\n? const String_["<<Idx<<"]: Out of range, [0] used\n";
	Idx=0;
    }
    return(S[Idx]);
}

char& String_::operator[](int Idx)
{
    if (Idx<0 || Idx>=max_len())
    {
	cerr<<"\n? String_["<<Idx<<"]: Out of range, [0] used\n";
	Idx=0;
    }
    return(S[Idx]);
}
// END of []

/* strchr(),strrchr(): returns the index of the first/last occurrence
 * of char C in the calling object or -1 if not found. Note that this
 * behaviour is different from the C functions which return pointers.
 */
int String_::strchr(int C) const
{
    char *Pos=::strchr(S, C);
    return((Pos==NULL)? -1: Pos-S);
}

int String_::strrchr(int C) const
{
    char *Pos=::strrchr(S, C);
    return((Pos==NULL)? -1: Pos-S);
}
// END of str[r]chr()

/* strstr(): returns the index of the first occurrence of Str in the
 * calling object or -1 if it was not there. Note the difference between
 * this and the ANSI C behaviour which returns a pointer.
 */
int String_::strstr(const String_& Str) const
{
    char *Pos=::strstr(S, Str.S);
    return((Pos==NULL)? -1: Pos-S);
}
// END of strstr()

/* max_len(L): sets the max. length of the string to L and returns the old length.
 * If the new length is longer than the old, then the extra new chars
 * will be set to '\0', if shorter, then the overhang will be truncated, 
 * if L is the same as the old length then nothing happens.
 */
unsigned int String_::max_len(unsigned int L)
{
    unsigned int Olen=Len;
    if (++L==Olen) return(Olen);	// do nothing
    
    char *Snew=new char [L];
    memcpy(Snew, S, (L>Olen)? Olen: L);
    if (Olen>L) Snew[L-1]='\0'; else memset(Snew+Olen, 0, L-Olen);
    delete [] S; S=Snew; Len=L;
    return(Olen-1);
}
// END of max_len(L)

// ---- Concatenation ----

/* S+=Str: appends Str to S and returns S.
 * S+Str: appends Str to S and returns the resulting object.
 */
String_& String_::operator+=(const String_& Str)
{
    if (!Str) return(*this);	// Str was empty, do nothing
    
    unsigned int Ls=strlen(Str.S);
    char *Snew=new char [Len=strlen(S)+Ls+1];
    strcpy(Snew, S); strcat(Snew, Str.S);
    delete [] S; S=Snew; 
    return(*this);
}
// END of +=

String_ String_::operator+(const String_& Str) const
{
    String_ Stemp(*this);
    Stemp+=Str;
    return(Stemp);
}
// END of +

// ---- Input/Output ----

/* In>>Str: Skips whitespaces and reads non-whitespace chars into
 * Str until it gets filled up or a whitespace is encountered. Will NOT
 * adjust the length of Str dynamically! 
 */
istream& operator>>(istream& In, String_& Str)
{
    if (!In) return(In);    // error state was set, do nothing
    In>>setw(Str.max_len()+1)>>(Str.S);	// do the read
    return(In);
}
// END of >>

/* Out<<Str: writes the string Str to Out. */
ostream& operator<<(ostream& Out, const String_& Str)
{
    Out<<Str.S;
    return(Out);
}
// END of <<

// ==== END OF METHODS String.c++ ====
