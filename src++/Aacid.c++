// ==== PROJECT DRAGON: METHODS Aacid.c++ ====

/* Class for storing amino acid atomic coordinates. */

// SGI C++ 4.0, IRIX 5.3, 4. June 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

// ---- MODULE HEADER ----

#include "Aacid.h"

// ==== METHODS ====

// ---- Static init ----

const char* Aacid_::Aas="ACDEFGHIKLMNPQRSTVWY";

// ---- Access ----

/* find(): locates an atom Name in the calling object.
 * Returns -1 if not found, the "absolute" (unmasked)
 * index otherwise.
 */
inline
int Aacid_::find(const char *Name) const
{
    if (Name==NULL || strlen(Name)>4) return(-1); // can't be an atom name
    int i;
    for (i=0; i<len(); i++)
	if (!strncmp(Atnames+4*i, Name, 4)) return(i);
    return(-1);
}
// END of find()

/* active(): invoked with an atom Name, returns true if
 * Name was found in the calling object and was active, 
 * false otherwise. 
 * Sets the activation status of Name to Flag if Name
 * was found and returns true: does nothing and
 * returns false if Name was not found.
 */
bool Aacid_::active(const char *Name) const
{
    int I=find(Name);
    return((I>=0)? Points_::active(I): false);
}
bool Aacid_::active(const char *Name, bool Flag)
{
    int I=find(Name);
    return((I>=0)? Points_::active(I, Flag), true: false);
}
// END of active()

/* atom(): returns (const) ptr to the coords of an atom called Name or
 * NULL if Name is not found or inactive.
 */
const Vector_* Aacid_::atom(const char *Name) const
{
    int I=find(Name);
    return((I>=0)? Data+I: NULL);
}
Vector_* Aacid_::atom(const char *Name)
{
    int I=find(Name);
    return((I>=0)? Data+I: NULL);
}
// END of atom()

// ---- Main and side chain masking ----

/* main_chain(): without an argument, returns true if all four
 * main-chain atoms (N, CA, C, O) are active, false otherwise.
 * With a bool argument Flag, switches the main chain on (true) or off (false).
 * Returns previous status.
 */
bool Aacid_::main_chain() const
{
    return(bool(active("N") && active("CA") && active("C") && active("O")));
}
bool Aacid_::main_chain(bool Flag)
{
    bool Oldstatus=main_chain();
    active("N", Flag); active("CA", Flag); 
    active("C", Flag); active("O", Flag);
    return(Oldstatus);
}
// END of main_chain()

/* side_chain(): without an argument, returns true if all
 * side-chain atoms are active, false otherwise (always false for Gly).
 * With a bool argument Flag, switches the complete side chain on
 * (true) or off (false), has no effect on Gly. Returns previous status.
 */
bool Aacid_::side_chain() const
{
    if (res_id()=='G') return(false);
    for (register unsigned int i=4; i<len(); i++)
	if (NULL==name(i)) return(false);
    return(true);
}

bool Aacid_::side_chain(bool Flag)
{
    if (res_id()=='G') return(false);
    bool Oldstatus=side_chain();
    for (register unsigned int i=4; i<len(); i++)
	active(name(i), Flag);
    return(Oldstatus);
}
// END of side_chain()

// ---- Auxiliaries ----

/* check_aa(): checks if the char Aa denotes one of the standard
 * amino acids. Lowercase converted to uppercase, silly chars
 * are converted to 'G'. Static
 */
inline
char Aacid_::check_aa(char Aa)
{
    Aa=toupper(Aa);
    return(strchr(Aas, Aa)? Aa: 'G');
}
// END of check_aa()

/* res_id(): without an arg, returns the current residue type.
 * With a char argument, sets the type of the calling object
 * to the corresponding amino acid or Gly if the char was not
 * a valid 1-letter amino acid code. Returns old type.
 */
char Aacid_::res_id(char Aa)
{
    char Oldaa=Restype; Restype=check_aa(Aa);
    setup(Restype); return(Oldaa);
}
// END of res_id()

/* alloc_array(): readjusts storage size to Size. Deletes
 * old arrays and allocates new ones. The tricky bit is the
 * Atnames array which contains the atom names as a contiguous string
 * made up of the concatenation of 4-letter "pdbprot" left-justified
 * atom IDs. All points are switched on.
 */
inline
void Aacid_::alloc_array(unsigned int Size)
{
    if (Size!=len())	// did not fit, realloc needed
    {
	delete [] Atnames;
	Atnames=new char [4*Size+1];
	len(Size);
    }
    Atnames[0]='\0'; mask(true);
}
// END of alloc_array()

/* copy_data(): copies the contents of A into the calling object.
 * Used by the copy ctor and the assignment. 
 */
inline
void Aacid_::copy_data(const Aacid_& A)
{
    setup(A.res_id());	// do the realloc and the atom descriptions
    Points_::operator=(A);	// assign base class
}
// END of copy_data()

Aacid_& Aacid_::operator=(const Aacid_& A)
{
    if (this==&A) return(*this);	// x=x
    copy_data(A);
    return(*this);
}
// =

/* setup(): deletes previous data structure and sets up the calling
 * object so that it can store all atoms for amino acid Aa. Aa is
 * supposed to have been checked for validity (default G). The new data array
 * will contain undefined values.
 */
void Aacid_::setup(char Aa)
{
    // set up new structure
    switch (Aa)
    {
	case 'A':
	alloc_array(5); memcpy(Atnames+16, "CB\0\0", 4);
	break;
	
	case 'C':
	case 'S':
	alloc_array(6);
	memcpy(Atnames+16, (Aa=='C')? "CB\0\0SG\0\0": "CB\0\0OG\0\0", 8);
	break;
	
	case 'D':
	case 'N':
	alloc_array(8); 
	memcpy(Atnames+16, (Aa=='D')? "CB\0\0CG\0\0OD1\0OD2\0": "CB\0\0CG\0\0OD1\0ND2\0", 16);
	break;
	
	case 'E':
	case 'Q':
	alloc_array(9); 
	memcpy(Atnames+16, (Aa=='E')? "CB\0\0CG\0\0CD\0\0OE1\0OE2\0": "CB\0\0CG\0\0CD\0\0OE1\0NE2\0", 20);
	break;
	
	case 'F':
	case 'Y':
	if (Aa=='F')
	{
	    alloc_array(11);
	    memcpy(Atnames+16, "CB\0\0CG\0\0CD1\0CD2\0CE1\0CE2\0CZ\0\0", 28);
	}
	else
	{
	    alloc_array(12);
	    memcpy(Atnames+16, "CB\0\0CG\0\0CD1\0CD2\0CE1\0CE2\0CZ\0\0OH\0\0", 32);
	}
	break;
	
	case 'H':
	alloc_array(10); memcpy(Atnames+16, "CB\0\0CG\0\0ND1\0CD2\0CE1\0NE2\0", 24);
	break;
	
	case 'I':
	alloc_array(8); memcpy(Atnames+16, "CB\0\0CG1\0CG2\0CD1\0", 16);
	break;
	
	case 'K':
	alloc_array(9); memcpy(Atnames+16, "CB\0\0CG\0\0CD\0\0CE\0\0NZ\0\0", 20);
	break;
	
	case 'L':
	alloc_array(8); memcpy(Atnames+16, "CB\0\0CG\0\0CD1\0CD2\0", 16);
	break;
	
	case 'M':
	alloc_array(8); memcpy(Atnames+16, "CB\0\0CG\0\0SD\0\0CE\0\0", 16);
	break;
	
	case 'P':
	alloc_array(7); memcpy(Atnames+16, "CB\0\0CG\0\0CD\0\0", 12);
	break;
	
	case 'R':
	alloc_array(11); memcpy(Atnames+16, "CB\0\0CG\0\0CD\0\0NE\0\0CZ\0\0NH1\0NH2\0", 28);
	break;
	
	case 'T':
	case 'V':
	alloc_array(7);
	memcpy(Atnames+16, (Aa=='T')? "CB\0\0OG1\0CG2\0": "CB\0\0CG1\0CG2\0", 12);
	break;
	
	case 'W':
	alloc_array(14);
	memcpy(Atnames+16, "CB\0\0CG\0\0CD1\0CD2\0NE1\0CE2\0CE3\0CZ2\0CZ3\0CH2\0", 40);
	break;
	
	case 'G':   // any unknown species has main-chain only
	default:
	alloc_array(4);
	break;
    }	    // switch Aa
    
    // construct mainchain bit
    memcpy(Atnames, "N\0\0\0CA\0\0C\0\0\0O\0\0\0", 16);
    Restype=Aa;
}
// END of setup()

// ---- Output ----

ostream& operator<<(ostream& Out, const Aacid_& A)
{
    Out<<"Residue=\'"<<A.res_id()<<"\'\n";
    const char *Name=NULL;
    for (int i=0; i<A.atom_no(); i++)
    {
	Name=A.name(i); Out<<Name;
	if (A.active(Name))
	    Out<<" +\n";    // ##### <<(*(A.atom(Name)));
	else
	    Out<<" -\n";
    }
    return(Out);
}
// END of <<

// ==== END OF METHODS Aacid.c++ ====
