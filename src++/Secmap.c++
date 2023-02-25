// ==== CLASS METHODS Secmap.c++ ====

/* Classes for mapping secondary structure information
 * from known structures onto target sequences via
 * a multiple alignment.
 */

// SGI C++, 27-Jun-1998. Andras Aszodi

// ---- MODULE HEADER ----

#include "Secmap.h"

// ==== METHODS ====

// ---- Smap_ METHODS ----

// set_nonbeta(): set the calling object to the secondary structure
// type specified by Stype. Prints a warning and does no assignment
// if Stype==BETA: use set_beta() for beta assignments.
// Return value: 1 for success, 0 on failure.
int Smap_::set_nonbeta(Sectype_ Stype)
{
  if (Stype==BETA)
  {
    cerr<<"\n? Smap_::set_nonbeta(): Use \'set_beta()\' for beta assignments"<<endl;
    return(0);
  }

  Sectype=Stype;
  reset_nonbeta();
  return(1);
}
// END of set_nonbeta()

// set_beta(): set the calling object to beta structure assignment.
// Specify the partner residues via Pn1,Pn2 (0 if there is no partner,
// ie. edge strand: however, Pn1 and Pn2 cannot be 0 at the same time)
// and the direction of the partner strands in A1,A2 (true for antiparallel,
// false for parallel: if the corresponding Pn is 0 then A is ignored).
// ID is the sheet ID character from DSSP.
// Return value: 0 on error, 1 if OK.
int Smap_::set_beta(int Pn1, int Pn2, bool A1, bool A2, char ID)
{
  if (!Pn1 && !Pn2)
  {
    cerr<<"\n? Smap_::set_beta(): Both beta partners are 0\n";
    return(0);
  }
  if (Pn1<0)
  {
    cerr<<"\n? Smap_::set_beta(): Partner 1 is "<<Pn1<<"<0: abs val used\n";
    Pn1=-Pn1;
  }
  if (Pn2<0)
  {
    cerr<<"\n? Smap_::set_beta(): Partner 2 is "<<Pn2<<"<0: abs val used\n";
    Pn2=-Pn2;
  }
  Sectype=BETA;
  Partner1=Pn1; Anti1=A1;
  Partner2=Pn2; Anti2=A2;
  Sheetid=ID;
  return(1);
}
// END of set_beta()

// get_betapartner(): when invoked for beta-containing objects,
// then the partner residue number is returned in Pn, and the
// partner direction in A (true for antiparallel, false for parallel).
// If Partnerno <=1, then Partner 1 data are returned, otherwise Partner 2.
// Return value: 0 on error (if the calling object is non-beta), 1 if OK.
int Smap_::get_betapartner(int Partnerno, int& Pn, bool& A) const
{
  if (Sectype!=BETA)
  {
    cerr<<"Smap_::get_betapartner(): Sectype is not BETA\n";
    return(0);
  }
  if (Partnerno<=1)
  {
    Pn=Partner1; A=Anti1;
  }
  else
  {
    Pn=Partner2; A=Anti2;
  }
  return(1);
}
// END of get_betapartner()

// Equivalence: two Smap_ objects are equal either if they
// both have the same non-BETA sectype or when they are both
// BETA and the partner info matches as well.
bool Smap_::operator==(const Smap_& S) const
{
  if (sec_type()!=S.sec_type()) return(false);
  else
  {
    if (sec_type()==BETA)
      return(Partner1==S.Partner1 && Anti1==S.Anti1 &&
	     Partner2==S.Partner2 && Anti2==S.Anti2);
    else return(true);
  }
}
// END of ==

// reset_nonbeta(): resets the beta-related members of the calling object
// to zeros and falses etc. Prints a warning and does no assignment
// if the calling object has a BETA secstr status.
// Return value: 1 on success, 0 on failure.
int Smap_::reset_nonbeta()
{
  if (Sectype==BETA)
  {
    cerr<<"\n? Smap_::reset_nonbeta(): Invoke me for non-beta objects only"<<endl;
    return(0);
  }
  
  Anti1=Anti2=false;
  Partner1=Partner2=0;
  return(1);
}
// END of reset_nonbeta()

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
ostream& operator<<(ostream& Out, const Smap_& Smap)
{
  if (Smap.sec_type()==Smap_::BETA)
  {
    int Pn;
    bool A;

    Smap.get_betapartner(1, Pn, A);
    if (!Pn) Out<<"----";
    else Out<<setw(3)<<Pn<<(A? 'a': 'p');
    Out<<Smap.Sheetid;
    Smap.get_betapartner(2, Pn, A);
    if (!Pn) Out<<"----";
    else Out<<setw(3)<<Pn<<(A? 'a': 'p');
  }
  else
  {
    Out<<"    ";
    switch(Smap.sec_type())
    {
      case Smap_::GAP: Out<<'-'; break;
      case Smap_::HELIX_310: Out<<'3'; break;
      case Smap_::HELIX_AL: Out<<'h'; break;
      case Smap_::HELIX_PI: Out<<'p'; break;
      case Smap_::OTHER: default: Out<<' '; break;
    }
    Out<<"    ";
  }
  return(Out);
}
// END of <<

// ---- Secmap_ METHODS ----

// Assignment
Secmap_& Secmap_::operator=(const Secmap_& Sm)
{
  if (this==&Sm) return(*this);  // x=x
  if (Sm.Secsno>Secsno)  // must grow
  {
    delete [] Secs;
    Secs=new Smap_ [Sm.Secsno];
  }
  Aa=Sm.Aa; Resno=Sm.Resno; Secsno=Sm.Secsno; 
  Cons=Sm.Cons; Secons=Sm.Secons;
  for (int i=0; i<Secsno; i++) Secs[i]=Sm.Secs[i];
  return(*this);
}
// END of =

// []: returns a const reference to the Idx:th secstr mapping.
// If Idx<0, then 0 is used, if Idx>=Secsno, then Secsno-1.
// Note that there is no non-const [] access.
const Smap_& Secmap_::operator[](int Idx) const
{
  if (Idx<0)
  {
    cerr<<"\n? Secmap_["<<Idx<<"]: returning [0]\n";
    Idx=0;
  }
  if (Idx>=Secsno)
  {
    cerr<<"\n? Secmap_["<<Idx<<"]: returning ["<<(Secsno-1)<<"]\n";
    Idx=Secsno-1;
  }
  return(Secs[Idx]);
}
// END of []

// set_aa(): sets the amino acid code to Ax,
// residue number to Rn, number of scaffold structures to Sn,
// and the mapping to a consistent OTHER.
void Secmap_::set_aa(char Ax, int Rn, int Sn)
{
  Aa=Ax; Resno=Rn;
  Cons=true; Secons=Smap_();
  if (Sn>Secsno)  // must grow
  {
	delete [] Secs;
	Secs=new Smap_ [Sn];
  }
    Secsno=Sn;
    for (int i=0; i<Secsno; i++) Secs[i]=Secons;
}
// END of set_aa()

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
int Secmap_::set_struct(int Idx, const Smap_& Smap)
{
  if (Idx<0 || Idx>=Secsno)
  {
    cerr<<"\n? Secmap_::sec_struct(): Index "<<Idx<<" outside range [0.."
        <<(Secsno-1)<<"]\n";
    return(0);
  }

  // check mapping consistency
  if (Cons)
  {
    // was consistent so far; this can be destroyed
    // only if the new mapping and the old consistent
    // mapping are both non-OTHER and they are not
    // equivalent
    if (Secons.sec_type()!=Smap_::OTHER && Smap.sec_type()!=Smap_::OTHER
	  && Secons!=Smap)
      Cons=false;
  }
  else
  {
    // was not consistent; this status remains if the new
    // mapping in Smap is the same as the old.
    // If it is different from the old, then a full check
    // against the others is necessary
    if (Smap==Secs[Idx]) Cons=false;
    else
    {
      bool C=true;
      Secons=Smap_();  // set to OTHER
      for (int i=0; C && i<Secsno; i++)
      {
	if (i==Idx) continue;  // no check against old value
	if (Secons.sec_type()==Smap_::OTHER && Secs[i].sec_type()!=Smap_::OTHER)
	{
	  Secons=Secs[i];  // overrides OTHER (first helix or beta)
	  continue;        // but still consistent
	}
	C=(Secons==Secs[i]);  // check consistency here
      }
      Cons=C && (Secons==Smap);
    }
  }
  Secs[Idx]=Smap;

    // non-OTHER mappings override the consistent OTHER
    if (Cons && Secons.sec_type()==Smap_::OTHER && Smap.sec_type()!=Smap_::OTHER)
	Secons=Smap;
  return(1);
}
// END of set_struct()

// <<: lists the calling object to Out. The format is the following line:-
// "X [rno] consmap map1 map2 ... mapN\n"
// where X is the 1-letter amino acid code, 'rno' is the residue number,
// consmap is the consistent secondary structure mapping or the string
// "non-cons!" if the mapping is not consistent. 'map1' ... 'mapN' are
// the individual mappings from the N scaffolds: for the format see
// the << operator for the Smap_ class.
ostream& operator<<(ostream& Out, const Secmap_& S)
{
  Out<<S.aa();
  if (S.aa()=='-' || !S.Resno) Out<<" [---] ";
  else Out<<" ["<<setw(3)<<S.Resno<<"] ";
  if (S.cons()) Out<<S.Secons; else Out<<"non-cons!";
  for (int i=0; i<S.Secsno; i++) Out<<' '<<S.Secs[i];
  Out<<endl;
  return(Out);
}
// END of <<

// ==== END OF METHODS Secmap.c++ ====

