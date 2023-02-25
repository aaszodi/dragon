#ifndef SECMAP_HEADER
#define SECMAP_HEADER

// ==== CLASS HEADER Secmap.h ====

/* Classes for mapping secondary structure information
 * from known structures onto target sequences via
 * a multiple alignment.
 */

// SGI C++, 27-Jun-1998. Andras Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

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

// Class Smap_: contains the description of a
// secondary structure assignment for a given target amino acid
// which comes from a known structure in the alignment. 
class Smap_
{
    // enum for secstr type
    public:
    enum Sectype_ { GAP=0, HELIX_310, HELIX_AL, HELIX_PI, BETA, OTHER };
    
    // data
    private:
    Sectype_ Sectype;	// secondary structure type
    bool Anti1, Anti2;	// true if antipar beta partner, false if parallel
    int Partner1, Partner2; // beta partner AA no in TARGET, 0 if unpaired
    char Sheetid;   	// beta-sheet ID (always CAPS), " " if non-beta
    
    // methods
    public:
    
    // default ctor
    Smap_(): Sectype(OTHER) { reset_nonbeta(); }
    
    // access

    // set_nonbeta(): set the calling object to the secondary structure
    // type specified by Stype. Prints a warning and does no assignment
    // if Stype==BETA: use set_beta() for beta assignments.
    // Return value: 1 for success, 0 on failure.
    int set_nonbeta(Sectype_ Stype);
    
    // set_beta(): set the calling object to beta structure assignment.
    // Specify the partner residues via Pn1,Pn2 (0 if there is no partner,
    // ie. edge strand: however, Pn1 and Pn2 cannot be 0 at the same time)
    // and the direction of the partner strands in A1,A2 (true for antiparallel,
    // false for parallel: if the corresponding Pn is 0 then A is ignored).
    // ID is the sheet ID character from DSSP.
    // Return value: 0 on error, 1 if OK.
    int set_beta(int Pn1, int Pn2, bool A1, bool A2, char ID);
    
    Sectype_ sec_type() const { return(Sectype); }
    
    // get_betapartner(): when invoked for beta-containing objects,
    // then the partner residue number is returned in Pn, and the
    // partner direction in A (true for antiparallel, false for parallel).
    // If Partnerno <=1, then Partner 1 data are returned, otherwise Partner 2.
    // Return value: 0 on error (if the calling object is non-beta), 1 if OK.
    int get_betapartner(int Partnerno, int& Pn, bool& A) const;
    
    // Equivalence: two Smap_ objects are equal either if they
    // both have the same non-BETA sectype or when they are both
    // BETA and the partner info matches as well.
    bool operator==(const Smap_& S) const;
    bool operator!=(const Smap_& S) const { return(!(*this==S)); }

    // <<: prints the calling object. For beta assignments, the format is
    // 9 chars wide and looks like "xxxpBpyyy" where xxx and yyy are the
    // partner residue numbers (or '---' for no partner), p is the char 'p'
    // for parallel partners 'a' for antiparallel partners, or '-' if there
    // is no partner, 'B' is the beta-sheet ID. For the other assignments,
    // there is one character printed in the middle of the 9-char field
    // so that it aligns nicely with the sheet IDs. The symbols are:-
    // '-' for gap, '3' for 3/10 helix, 'h' for alpha-helix, 'p' for pi-helix,
    // ' ' for "other" (coil or unspecified etc).
    // There is no return char printed.
    friend ostream& operator<<(ostream& Out, const Smap_& Smap);
    
    protected:
    int reset_nonbeta();
};
// END OF CLASS Smap_

// Class Secmap_: contains a description of all secondary structure
// assignments obtained from the known scaffold structures
// for a given residue in the target sequence.
class Secmap_
{
  // data
  protected:

  char Aa;  // 1-letter amino acid code
  int Resno;  // residue number in target sequence
  Smap_ *Secs;  // array of secstr assignments from scaffolds
  int Secsno;   // number of scaffolds (ie. secstr assignments)
  bool Cons;  // true if the assignments are consistent (same from all)
  Smap_ Secons; // if Cons==true, this holds the consistent secstr assignment

  // methods
  public:

  // default ctor: init to hold N assignments for amino acid Ax, resno Rn.
  // No validity checks on params
  Secmap_(char Ax='X', int Rn=0, int N=0):
    Aa(Ax), Resno(Rn), Secs(NULL), Secsno(N), Cons(true), Secons()
  {
    if (Secsno) Secs=new Smap_ [Secsno];
  }

  // copy ctor
  Secmap_(const Secmap_& Sm):
    Aa(Sm.Aa), Resno(Sm.Resno), Secs(NULL), Secsno(Sm.Secsno), 
    Cons(Sm.Cons), Secons(Sm.Secons)
  {
    if (Secsno)
    {
      Secs=new Smap_ [Secsno];
      for (int i=0; i<Secsno; i++) Secs[i]=Sm.Secs[i];
    }
  }

  // dtor
  ~Secmap_() { delete [] Secs; }

  // Assignment
  Secmap_& operator=(const Secmap_& Sm);
  
  // simple access
  char aa() const { return(Aa); }
  int resno() const { return(Resno); }
  bool cons() const { return(Cons); }
  
    // []: returns a const reference to the Idx:th secstr mapping.
    // If Idx<0, then 0 is used, if Idx>=Secsno, then Secsno-1.
    // Note that there is no non-const [] access.
    const Smap_& operator[](int Idx) const;

    // set_aa(): sets the amino acid code to Ax,
    // residue number to Rn, number of scaffold structures to Sn,
    // and the mapping to a consistent OTHER.
    void set_aa(char Ax, int Rn, int Sn);

  // set_struct(): sets the Idx:th scaffold assignment to
  // what is contained in Smap. No assignment is done if
  // Idx is outside the allowed index range. If OK, then
  // the secstr mapping is stored and is compared to the
  // mappings already in the calling object to check the consensus.
  // This means that a consensus secstr mapping may have
  // only OTHER and/or one of the remaining secstr types,
  // and if the non-OTHER type is BETA, then all partner information
  // must match as well.
  // Return value: 0 on error, 1 if OK.
  int set_struct(int Idx, const Smap_& Smap);

  // cons_struct(): returns true and the secstr mapping
  // in Smap if the mapping in the calling object was consistent,
  // false otherwise (in which case Smap is not modified).
  bool cons_struct(Smap_& Smap) const
  {
    if (Cons) Smap=Secons;
    return(Cons);
  }

    // <<: lists the calling object to Out. The format is the following line:-
    // "X [rno] consmap map1 map2 ... mapN\n"
    // where X is the 1-letter amino acid code, 'rno' is the residue number,
    // consmap is the consistent secondary structure mapping or the string
    // "non-cons!" if the mapping is not consistent. 'map1' ... 'mapN' are
    // the individual mappings from the N scaffolds: for the format see
    // the << operator for the Smap_ class.
    friend ostream& operator<<(ostream& Out, const Secmap_& S);

};
// END OF CLASS Secmap_

// ==== END OF HEADER Secmap.h ====

#endif  /* SECMAP_HEADER */

