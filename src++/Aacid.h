#ifndef AACID_CLASS
#define AACID_CLASS

// ==== PROJECT DRAGON: HEADER Aacid.h ====

/* Class for storing amino acid atomic coordinates. */

// SGI C++ 4.0, IRIX 5.3, 4. June 1996. Andris Aszodi

// ---- UTILITY HEADERS ----

#include "Points.h"

// ==== CLASSES ====

/* Aacid_: a class for storing amino acid descriptions
 * and data. An atom may be identified by its name such as "CG2".
 * Only atom names consistent with the amino acid type
 * in the object can be used. Derived from Points_.
 */
class Aacid_: public Points_
{
    // data
    protected:
    
    static const char* Aas; // legal 1-letter AA codes
    char Restype;   // 1-letter AA code
    char *Atnames;  // atom name array (contig. 4-letter IDs, \0-termination at end)
    
    // methods
    public:
    
	// ctors, dtor, assignment
    /* Inits to hold amino acid Aa (default 'G'). Data left uninitialised. */
    Aacid_(char Aa='G'): Points_(0), 
	Restype(Aa), Atnames(NULL) { setup(Aa); }
    
    /* Inits to hold the same data as A. */
    Aacid_(const Aacid_& A): Points_(A), 
	Restype(A.res_id())
    {
	Atnames=new char [4*len()+1];
	memcpy(Atnames, A.Atnames, 4*len()+1);
    }
    
    virtual ~Aacid_() { delete [] Atnames; }
    
    Aacid_& operator=(const Aacid_& A);
    
	// access
    /* res_id(): without an arg, returns the current residue type.
     * With a char argument, sets the type of the calling object
     * to the corresponding amino acid or Gly if the char was not
     * a valid 1-letter amino acid code. Returns old type.
     */
    char res_id() const { return(Restype); }
    char res_id(char Aa);
    
    /* active(): invoked with an atom Name, returns true if
     * Name was found in the calling object and was active, 
     * false otherwise. 
     * Sets the activation status of Name to Flag if Name
     * was found and returns true: does nothing and
     * returns false if Name was not found.
     */
    bool active(const char *Name) const;
    bool active(const char *Name, bool Flag);
    
    /* atom(): returns (const) ptr to the coords of an atom called Name or
     * NULL if Name is not found or inactive.
     */
    const Vector_* atom(const char *Name) const;
    Vector_* atom(const char *Name);
    
    /* name(): returns a const ptr to the Index:th atom name or 
     * NULL if Index is out of range. Activity status is ignored.
     */
    const char *name(unsigned int Index) const { return((Index<len())? Atnames+4*Index: NULL); }
    
    /* atom_no(): returns the number of atoms in the object. Just syntactic sugar */
    unsigned int atom_no() const { return(len()); }
    
	// main and side chain masking
    /* main_chain(): without an argument, returns true if all four
     * main-chain atoms (N, CA, C, O) are active, false otherwise.
     * With a bool argument Flag, switches the main chain on (true) or off (false).
     * Returns previous status.
     */
    bool main_chain() const;
    bool main_chain(bool Flag);

    /* side_chain(): without an argument, returns true if all
     * side-chain atoms are active, false otherwise (always false for Gly).
     * With a bool argument Flag, switches the complete side chain on
     * (true) or off (false), has no effect on Gly. Returns previous status.
     */
    bool side_chain() const;
    bool side_chain(bool Flag);

    protected:
    
    int find(const char *Name) const;
    static char check_aa(char Aa);
    void copy_data(const Aacid_& A);
    void setup(char Aa='G');
    void alloc_array(unsigned int Size);
    
};
// END OF CLASS Aacid_

ostream& operator<<(ostream& Out, const Aacid_& A);

// ==== END OF HEADER Aacid.h ====

#endif	/* AACID_CLASS */
