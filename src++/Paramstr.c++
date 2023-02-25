// ==== PROJECT DRAGON: METHODS Paramstr.c++ ====

/* Stores global string (filename) parameters. */

// SGI C++ 7.1, IRIX 6.2, 7. Mar. 1997. Andris Aszodi

// ---- MODULE HEADER ----

#include "Paramstr.h"

// ==== Paramstr_ METHODS ====

// ---- Constructors ----

/* Sets the maximal length of the string variables within to Size.
 * If Size==0 (default), then the Defval string's length is
 * used, if Size is less than the maximal length of Defval,
 * then it is adjusted.
 */
Paramstr_::Paramstr_(const char *Defval, unsigned int Size, 
	    const char *Nm, const char *Ds):
	Value(Defval), Default(Defval), Parambase_(Nm, Ds)
{
    if (Size)
    {
	if (Size<Default.max_len()) Value.max_len(Default.max_len());
	else
	{
	    Value.max_len(Size); Default.max_len(Size);
	}
    }
}

// ---- Access ----

/* Converts the calling object to a (const) char ptr.
 * If the value of the calling object begins with the string
 * "$DRAGON_DATA/", then this string will be substituted
 * with the value of the environment variable DRAGON_DATA
 * if defined, otherwise it is removed.
 */
Paramstr_::operator const char* () const
{
    const char *V=Value;
    static const int MACROLEN=13;
    static String_ Newval;  // to avoid nasty memory leaks
    
    if (!strncmp(V, "$DRAGON_DATA/", MACROLEN))
    {
	// prefix found
	char *Data=getenv("DRAGON_DATA");
	Newval=(Data==NULL)? ".": Data;
	Newval+=V+MACROLEN-1;
	return((const char*) Newval);
    }
    else return(V);
}
// END of const char* conversion

/* set_default(): sets the default string to Defval and the maximal
 * size to Size (default==0). If Size==0, then the length of Defval
 * is used as the maximal size (same as in the ctor). Sets the
 * Value string to Defval also.
 * Returns maximal size.
 */
unsigned int Paramstr_::set_default(const char *Defval, unsigned int Size)
{
    Value=Default=String_(Defval);
    if (Size)
    {
	if (Size<Default.max_len()) Size=Default.max_len();
	Value.max_len(Size); Default.max_len(Size);
    }
    else Size=Default.max_len(); Changed=true;
    return(Size);
}
// END of set_default()

// ---- Input/output ----

/* read_from(): reads the following line from an input stream In:
 * "NAME string"
 * where NAME should match the Name member of the object and "string"
 * should be a string not containing whitespaces. 
 * If the line was of the form "NAME " then the default string is
 * used as the value (which would have been "" otherwise).
 * If "string" is longer than the max. size in the object, then
 * the extra chars are left in In but this is OK.
 * Return values: -1 on error, 0 if NAME did not match, 1 if OK.
 */
int Paramstr_::read_from(istream& In)
{
    if (!In.good()) return(-1);
    
    String_ Namebuf(Name);
    Namebuf[0]='\0'; // make sure sizes match
    streampos Pos=In.tellg();	// save initial position
    
    In>>Namebuf;
    if (!In.good() || Namebuf!=Name)
    {
	In.seekg(Pos);	// name did not match, reset and return
	return(0);
    }
    
    In.flags(ios::skipws);
    In>>Value;
    if (!In)
    {
	cerr<<"\n? Paramstr_::read_from(): Cannot read string "<<Name<<endl;
	Value=Default;
	In.clear(ios::failbit|In.rdstate());
	return(-1);
    }
    if (!Value) Value=Default;	// "" is not allowed unless it is the default
    Changed=true;	// make a note that the value has changed 
    return(1);
}
// END of read_from()

/* Lists to Out. */
ostream& Paramstr_::write_to(ostream& Out, bool Comments) const
{
    if (Comments) Out<<"\n# "<<Descr<<" (default="<<Default<<")\n";
    Out<<Name<<" "<<Value<<"\n";
    return(Out);
}
// END of write_to()

// ==== END OF METHODS Paramstr.c++ ====
