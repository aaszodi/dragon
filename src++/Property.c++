// ==== PROJECT DRAGON: METHODS Property.c++ ====

/* Amino-acid property data storage. */

// SGI C++, IRIX 6.2, 14. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <strstream.h>
#include <fstream.h>
#include <iomanip.h>
#include <string.h>
#include <ctype.h>

// ---- CLASS HEADER ----

#include "Property.h"

// ---- DEFINITIONS ----

#define LINELEN 240

// ==== Property_ METHODS ====

// ---- Static data initialisation ----

// Amino acids in alphabetical order plus the GAP
const char Property_::GAP='-';
const char Property_::Aacodes[AANO+1]="ABCDEFGHIKLMNPQRSTVWXYZ-";

// Hydrophobicity values
const Proparray_ Property_::Hyphobdef={
    1.73, 0.02, 0.84, 0.03, 0.01, 1.48, 1.27, 0.06, 3.46, 
    0.03, 2.56, 0.86, 0.01, 0.18, 0.03, 0.00, 0.49, 0.59, 
    2.46, 0.74, 0.5, 0.59, 0.02, 0.0};

// Amino acid side chain volumes
const Proparray_  Property_::Volumedef={
    22.7, 50.2, 34.9, 46.5, 63.5, 91.1, 5.7, 74.7, 73.7, 
    79.5, 73.7, 74.6, 54.0, 45.3, 71.0, 100.4, 30.4, 47.4,
    56.7, 120.7, 50.0, 100.2, 67.2, 150.0};

// ---- Constructor ----

/* Inits so that the default values are taken from the external
 * location pointed to by Defptr. Prop is initialised to { 0.0 ... }
 * if Defptr==NULL (the default ;-)). Defptr is assumed to
 * point to something meaningful such as one of the static const
 * members above.
 */
Property_::Property_(const float *Defptr)
{
    if (Defptr==NULL) memset(Prop, 0, AANO*sizeof(float));
    else memcpy(Prop, Defptr, AANO*sizeof(float));
    Default=Defptr;
}

// ---- Access ----

/* reset(): resets to the default values in the array pointed to
 * by Default, or to {0.0 ...} if it was NULL (see ctor above).
 */
void Property_::reset()
{
    if (Default==NULL) memset(Prop, 0, AANO*sizeof(float));
    else memcpy(Prop, Default, AANO*sizeof(float));
}
// END of reset()

/* get_index: returns the position of the character Aa in Aacodes
 * or the index of 'X' times -1 if Aa was not found. 
 * Aa is converted to uppercase.
 */
int Property_::get_index(char Aa) const
{
    char *Pos;
    
    Aa=toupper(Aa);
    Pos=strchr(Aacodes, Aa);
    if (Pos==NULL)
    {
	Pos=strchr(Aacodes, 'X');
	return(Aacodes-Pos);
    }
    else return(Pos-Aacodes);
}
// END of get_index()

/* [X]: returns the property value for amino acid Aa. If Aa is
 * lowercase, then it will be converted to uppercase. Invalid
 * chars will be treated as 'X' (unknown) silently.
 */
float Property_::operator[](char Aa) const
{
    int Idx=get_index(Aa);
    return(Prop[(Idx<0)? -Idx: Idx]);
}
// END of []

/* avg_val(): returns the average value of the property for the
 * amino acids in the string Posstr (representing an alignment position).
 */
float Property_::avg_val(const char *Posstr) const
{
    unsigned int N;
    if (!(N=strlen(Posstr))) return(0.0);   // empty string
    
    register unsigned int i;
    float Avg=0.0;
    for (i=0; i<N; i++) Avg+=(*this)[Posstr[i]];
    Avg/=N;
    return(Avg);
}
// END of avg_val()

// ---- Input/output ----

/* read_file(): reads the property data file Fname. The format is:-
 * <char> <value> ...\n
 * where <char> must be a valid uppercase 1-letter amino acid code, 
 * and <value> is a floating-point number. Lines beginning with '#'
 * are comments and ignored. Non-parsable lines will be skipped
 * with a warning to cerr. 
 * NOTE: the format has changed from Version 3.x and the "title line"
 * is no longer supported. Put a '#' in front of the titles in the
 * old property data files to avoid error messages.
 * Return value: 0 on error, 1 if OK.
 */
int Property_::read_file(const char *Fname)
{
    ifstream Infile(Fname);
    if (!Infile)
    {
	cerr<<"\n? Property_::read_file("<<Fname<<"): Cannot open\n";
	return(0);
    }
    Infile>>(*this);
    Infile.close();
    return(Infile.good() || Infile.eof());
}
// END of read_file()

/* >>: read from stream In. See comments to read_file(). */
istream& operator>>(istream& In, Property_& Pr)
{
    if (!In)
    {
	cerr<<"\n? >>Property_: Cannot read from stream\n";
	return(In);
    }
    
    char Line[LINELEN];
    istrstream Inline(Line, LINELEN-1);
    int Idx, Lineno;
    char Aa;
    float Val;
    
    for (Lineno=0; In.good(); Lineno++)
    {
	In.getline(Line, LINELEN, '\n');
	if (Line[0]=='#') continue;
	
	Inline.seekg(ios::beg);
	Inline>>Aa>>Val;
	if (!Inline)
	{
	    cerr<<"\n? >>Property_: Cannot read line " <<(Lineno+1)<<endl;
	    continue;
	}
	
	Idx=Pr.get_index(Aa);
	if (Idx<0)
	{
	    cerr<<"\n? >>Property_: \'"<<Aa<<"\' unknown, skipped\n";
	    continue;
	}
	Pr.Prop[Idx]=Val;
    }
    return(In);
}
// END of >>

/* Output: list to stream Out neatly. */
ostream& operator<<(ostream& Out, const Property_& Pr)
{
    Out<<"AA    Value\n-----------\n";
    for (unsigned int i=0; i<AANO; i++)
	Out<<Pr.Aacodes[i]<<"    "<<Pr.Prop[i]<<endl;
    Out<<"-----------\n";
    return(Out);
}
// END of <<

// ==== END OF METHODS Property.c++
