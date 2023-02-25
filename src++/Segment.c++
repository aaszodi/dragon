// ==== PROJECT DRAGON: METHODS Segment.c++ ====

/* Classes for handling a model polypeptide as non-overlapping
 * segments (of secondary structures). These classes define
 * the chain topology only; H-bond related stuff is in Secstr.c++.
 */

// SGI C++, IRIX 6.2, 9. Aug. 1996. (C) Andras Aszodi

// ---- HEADER ----

#include "Segment.h"

// ==== Linsegm_ METHODS ====

// ---- Access ----

/* beg(), end(): without args, the current beginning and end are returned.
 * With an uint arg, the beginning or the end are set. Zero-length segments
 * and Beg>End cases are disallowed and no action is taken if these settings
 * are attempted. The previous limits are returned.
 */
inline unsigned int Linsegm_::beg(unsigned int Newbeg)
{
    unsigned int Oldbeg=Beg;
    if (Newbeg!=Beg && Newbeg<=End)
    {
	Beg=Newbeg; Changed=true;
    }
    return(Oldbeg);
}
// END of beg()

inline unsigned int Linsegm_::end(unsigned int Newend)
{
    unsigned int Oldend=End;
    if (Newend!=End && Newend>=Beg)
    {
	End=Newend; Changed=true;
    }
    return(Oldend);
}
// END of end() ;-)

/* limits(): sets the beginning to Newbeg and the end to Newend
 * simultaneously. (The methods beg(B) and end(E) above have built-in
 * checks for disabling specifications such as Newbeg>=End and
 * therefore are not suitable for complete modifications.)
 * If Newbeg>Newend then the values are swapped.
 */
void Linsegm_::limits(unsigned int Newbeg, unsigned int Newend)
{
    if (Newbeg>Newend)
    {
	unsigned int T=Newbeg; Newbeg=Newend; Newend=T;
    }
    if (Newbeg!=Beg) { Beg=Newbeg; Changed=true; }
    if (Newend!=End) { End=Newend; Changed=true; }
}
// END of limits()

// ---- Membership ----

/* mask(): constructs and returns a bitmap in which all positions
 * corresponding to the members of the segment are set true, the others
 * false. Since Linsegm_ objects know nothing about the model chain length, 
 * this should be supplied in Rno. If Rno==0 (the default) then the
 * mask will be end()+1 positions long which can be set to the correct
 * chain length later. If Rno>0 but is less than enough for containing the
 * segment then a warning is printed and an end()+1 -long mask is created.
 */
Bits_ Linsegm_::mask(unsigned int Rno) const
{
    if (Rno && End+1>Rno)
    {
	cerr<<"? Linsegm_::mask(): Chain length "<<Rno<<" too short for this segment\n";
	Rno=End+1;
    }
    if (!Rno) Rno=End+1;    // Rno "unknown" (default)
    Bits_ Mask(Rno);    // set mask to proper length, all bits false by default
    
    for (unsigned int i=Beg; i<=End; i++) Mask.set_bit(i);  // flick members on
    return(Mask);
}
// END of mask()

// ==== Strand_ METHODS ====

/* sense(): w/o arguments, returns the Sense value.
 * With an int argument S, sets the sense to ANTI if S<0, 
 * to NONE if S==0, PAR if S>0 and returns the old sense value.
 */
inline Strand_::Sense_ Strand_::sense(int S)
{
    Sense_ Oldsens=Sense;
    Sense=(S>0)? PAR: ((S<0)? ANTI: NONE);
    return(Oldsens);
}
// END of sense()

/* phase(): w/o arguments, returns the phase.
 * With an int argument P, sets the phase to P and returns the
 * old phase.
 */
inline int Strand_::phase(int P)
{
    int Oldphase=Phase; Phase=P; return(Oldphase);
}
// END of phase()

// ==== Sheet_ METHODS ====

// ---- Constructors ----

/* Inits with a single strand. Note that 1-strand sheets do not exist
 * so use this just as a starting point for building sheets.
 * Sets the first strand's phase to 0 and sense to NONE.
 */
Sheet_::Sheet_(const Strand_& Str1) : Strands(1)
{ 
    Strands[0]=Str1;
    Strands[0].sense(Strand_::NONE); 
    Strands[0].phase(0);
}

// ---- Access ----

/* [Sno]: returns the Sno-th strand of the sheet. If Sno is
 * invalid, a warning is printed and the 0-th strand is returned.
 * Note that the returned value cannot be an lvalue.
 */
const Strand_& Sheet_::operator[](unsigned int Sno) const
{
    if (Sno>=Strands.len())
    {
	cerr<<"? Sheet_::operator[]: Index "<<Sno<<" out of range, 0 used\n";
	Sno=0;
    }
    return(Strands[Sno]);
}
// END of operator []

// ---- Membership ----

/* member(): returns true if Resno is within the sheet, false otherwise. */
bool Sheet_::member(unsigned int Resno) const
{
    for (unsigned int i=0; i<Strands.len(); i++)
	if (Strands[i].member(Resno)) return(true);
    return(false);
}
// END of member()

/* mask(): returns a Bits_ object in which all the bits corresponding to
 * the members of the sheet are switched ON, the rest OFF. Since Sheet_
 * objects know nothing about the size of the chain, the chain length
 * Rno must be specified. If it is omitted, then a default Rno==0
 * is assumed which means that the Bits_ will have a length suitable
 * for accommodating all members. It can then be enlarged later.
 * If Rno>0 is specified but it is too small, then it is adjusted.
 */
Bits_ Sheet_::mask(unsigned int Rno) const
{
    // get the largest end-position from the strands
    unsigned int i, Maxend=0, Sno=Strands.len();
    for (i=0; i<Sno; i++) 
	if (Strands[i].end()>Maxend) Maxend=Strands[i].end();
    
    // compare w/ Rno
    if (Rno && Rno<=Maxend)
    {
	cerr<<"? Sheet_::mask(): Chain length "<<Rno<<" is too short, adjusted\n";
	Rno=Maxend+1;
    }
    if (!Rno) Rno=Maxend+1;	// "unknown" chain length
    
    // construct mask
    Bits_ Mask(Rno), Smask(Rno);
    for (i=0; i<Sno; i++)
    {
	Smask=Strands[i].mask(Rno);    // get the strand mask
	Mask|=Smask;	// OR them together
    }
    return(Mask);
}
// END of mask()

// ---- Sheet building ----

/* first_strand(): completely resets the sheet and adds Str1
 * as its first strand. The sense and phase information in Str1
 * is ignored: the phase is set to 0 and the sense to NONE inside.
 * If the sheet was initialised w/ a strand (see the relevant ctor)
 * then this method need not be used. 
 */
void Sheet_::first_strand(const Strand_& Str1)
{
    Strands.len(1); Strands[0]=Str1;
    Strands[0].sense(Strand_::NONE); 
    Strands[0].phase(0);
    Changed=true;
}
// END of first_strand()

/* add_strand(): adds a new strand Str to the calling object in
 * the orientation specified within, phased so that the Thisres-th
 * residue in the new strand matches the Otherres-th residue in
 * the last strand within. No action is taken if Str.sense()==NONE
 * or if Thisres and Otherres are not within their respective strands
 * or if the new strand overlaps w/ one of the old ones.
 * Return value: 0 on error, the new length if successful.
 */
int Sheet_::add_strand(const Strand_& Str, int Thisres, int Otherres)
{
    // checks
    if (Str.sense()==Strand_::NONE)
    {
	cerr<<"? Sheet_::add_strand(): Sense missing from new strand, not added\n";
	return(0);
    }
    if (!Str.member(Thisres))
    {
	cerr<<"? Sheet_::add_strand(): Residue "<<Thisres<<" is not in new strand\n";
	return(0);
    }
    int Sno=Strands.len();
    if (!Strands[Sno-1].member(Otherres))
    {
	cerr<<"? Sheet_::add_strand(): Residue "<<Otherres<<" is not in last strand\n";
	return(0);
    }
    for (unsigned int i=0; i<Sno; i++)
    {
	if (!(Strands[i].end()<Str.beg() || Str.end()<Strands[i].beg()))
	{
	    cerr<<"? Sheet_::add_strand(): New strand overlaps with strand "<<i<<endl;
	    return(0);
	}
    }
    
    Strands.len(Sno+1);  // make room for new strand
    Strands[Sno]=Str;    // store new
    
    // calc and store phase
    Strands[Sno].phase( (Str.sense()==Strand_::PAR)?
	Otherres-Thisres+Str.beg()-Strands[Sno-1].beg():
	Strands[Sno-1].end()-Otherres-Thisres+Str.beg() );
    
    Changed=true;
    return(Sno+1);  // success
}
// END of add_strand()

// ---- Strand positions ----

/* sense(): returns the relative sense of orientation between
 * the strands no. S1 and S2. Returns NONE if S1==S2 or if 
 * S1 or S2 are not valid strand numbers (with a warning).
 */
Strand_::Sense_ Sheet_::sense(unsigned int S1, unsigned int S2) const
{
    if (S1>=Strands.len())
    {
	cerr<<"? Sheet_::sense(): Invalid strand index "<<S1<<endl;
	return(Strand_::NONE);
    }
    if (S2>=Strands.len())
    {
	cerr<<"? Sheet_::sense(): Invalid strand index "<<S2<<endl;
	return(Strand_::NONE);
    }
    if (S1==S2) return(Strand_::NONE);
    
    unsigned int i;
    int Sn=1;
    
    if (S1>S2) { i=S1; S1=S2; S2=i; }	// swap so that S1<S2
    for (i=S2; i>S1; i--) Sn*=Strands[i].sense();
    return((Sn<0)? Strand_::ANTI : Strand_::PAR);
}
// END of sense()

/* strand_res(): returns the index of the strand which contains
 * the residue Resno. Returns -1 (NOT 0!!) if Resno is not in the sheet.
 * Functionally similar to the Linsegm_::member() method.
 */
int Sheet_::strand_res(unsigned int Resno) const
{
    int i;
    for (i=Strands.len()-1; i>=0 && !Strands[i].member(Resno); i--);
    return(i);
}
// END of strand_res()

/* offs_strd(): The overall sheet orientation
 * is defined by the first strand and its first res serves as the
 * origin of the 2D "sheet coord system". An offset Offs on strand Sno
 * is simply the difference of a position from the beginning of the strand, 
 * e.g. residue 25 has an offset +5 on a strand 20->30. It is not
 * required that Offs correspond to a "real" residue, i.e. any Offs value
 * is legal. However, Sno must be a legal strand no. This method
 * calculates the offset on Strand 0 of a position whose offset is
 * Offs on Strand Sno,  based on the interstrand phasing information.
 * Together with strand_res(), offs_strd() locates every residue within
 * the "sheet coord system". If Sno is invalid, a warning is printed 
 * and 0 is returned.
 */
int Sheet_::offs_strd(unsigned int Sno, int Offs) const
{
    if (Sno>=Strands.len())
    {
	cerr<<"? Sheet_::offs_strd(): Strandno="<<Sno<<" illegal, 0 returned\n";
	return(0);
    }
    
    // work our way back to Strand 0
    for ( ; Sno>0; Sno--)
    {
	if (Strands[Sno].sense()==Strand_::PAR)	    // parallel to prev
	    Offs+=Strands[Sno].phase();
	else					    // anti to prev
	    Offs=Strands[Sno-1].len()-1-(Offs+Strands[Sno].phase());
    }
    return(Offs);
}
// END of offs_strd()

// ==== END OF METHODS Segment.c++ ====
