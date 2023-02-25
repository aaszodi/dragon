#ifndef SEGMENT_CLASSES
#define SEGMENT_CLASSES

// ==== PROJECT DRAGON: HEADER Segment.h ====

/* Classes for handling a model polypeptide as non-overlapping
 * segments (of secondary structures). These classes define
 * the chain topology only; H-bond related stuff is in Secstr.c++.
 */

// SGI C++, IRIX 6.2, 9. Aug. 1996. (C) Andras Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "Bits.h"

// ==== CLASSES ====

/* The model chains are divided into segments which correspond to the
 * secondary structure layout in the current implementation. The classes
 * in this module represent the topology of the secstr segments, while
 * the classes in the "Secstr" module represent the geometry (ideal
 * distances and structure w/ chirality). 
 *
 * The inheritance graph of the whole family is as follows:-
 *
 *	         [Segmbase_]
 *                    |
 *	    +---------+----------+
 *	    :         :          :
 *	    V         V          V
 *	Linsegm_  [Sstrbase_]  Sheet_
 *         |           |          |
 *         +--------+  +--------+ |
 *         |        |  |        | |
 *	   V        V  V        V V
 *	Strand_     Helix_      Beta_
 * 
 * ( [] enclose abstract base classes, : denotes virtual ancestry)
 * 
 * Note that Sheet_ *contains* Strand_ objects to describe the strands
 * within the sheet.
 * This module implements the following classes:-
 *	Segmbase_: [ABC] for the whole lot
 *	Linsegm_ : a linear piece of contiguous residues on a chain
 *	Strand_: is a Linsegm_ representing beta-strands
 *	Sheet_ : a collection of Strand_ -s with par/anti and phase info.
 * The rest is implemented in the "Secstr" module.
 */

/* Class Segmbase_: the abstract base class for the whole hierarchy.
 * Declares the most basic "membership" fns
 * as pure virtuals. The only data member is a bool sentinel which
 * becomes true after non-const method calls ("dirty bit").
 */
class Segmbase_
{
    // data
    protected:
    
    bool Changed;	// will be set after size or layout changes
    
    // methods
    public:
	// default constructor
    Segmbase_(): Changed(true) {}
    
	// destructor
    virtual ~Segmbase_() {}
    
	// access
    /* strand_no(): returns the number of contiguous segments in the
     * calling object. (Meaningful for sheets only, 1 for everybody else.)
     */
    virtual unsigned int strand_no() const =0;
    
	// membership
    /* member(): will return true if Resno is in the object, false otherwise. */
    virtual bool member(unsigned int Resno) const =0;
    
    /* mask(): will construct and return a bitmap in which all positions
     * corresponding to the residues in the object are set true, the others
     * false. Since these objects know nothing about the model chain length, 
     * this should be supplied in Rno. If Rno==0 (the default) then the
     * mask will be long enough to contain all residues: an all-false tail
     * might be added later. This will happen also when Rno is not large enough.
     */
    virtual Bits_ mask(unsigned int Rno=0) const =0;

};
// END OF CLASS Segmbase_

/* Class Linsegm_ : implements a "linear segment" which is a stretch
 * of consecutive residues along the model chain. A linear segment can be
 * described by specifying its beginning and end. Our convention for
 * numbering residues will be C-style, beginning with residue 0. 
 * We also require that the beginning of a Linsegm_ be <= than its end.
 * Beg==End means a 1-long segment. A Linsegm_ can have its beginning
 * and end altered, and can be queried about its length or whether a
 * residue is a member of a segment. It can also supply a bitmask in which
 * the constituent residues are marked ON,  and the others OFF.
 * Derived from Segmbase_ .
 */
class Linsegm_ : public virtual Segmbase_
{
    // data
    protected:
    unsigned int Beg, End;
    
    // methods
    public:
    
	// constructor
    /* Inits a linear segment to begin at Start and end at Stop.
     * Disallows Start>Stop cases etc. The default segment starts and
     * ends at residue 0. No warnings are printed.
     * Note that there's no way to check the segment limits against the
     * length of the model chain :-( .
     */
    Linsegm_(unsigned int Start=0, unsigned int Stop=0)
    {
	if (Start>Stop) { Beg=Stop; End=Start; }
	else { Beg=Start; End=Stop; }
    }

	// access
    /* beg(), end(): without args, the current beginning and end are returned.
     * With an uint arg, the beginning or the end are set. Zero-length segments
     * and Beg>End cases are disallowed and no action is taken if these settings
     * are attempted. The previous limits are returned.
     */
    unsigned int beg() const { return(Beg); }
    unsigned int beg(unsigned int Newbeg);
    unsigned int end() const { return(End); }
    unsigned int end(unsigned int Newend);

    /* limits(): sets the beginning to Newbeg and the end to Newend
     * simultaneously. (The methods beg(B) and end(E) above have built-in
     * checks for disabling specifications such as Newbeg>=End and
     * therefore are not suitable for complete modifications.)
     * If Newbeg>Newend then the values are swapped.
     */
    void limits(unsigned int Newbeg, unsigned int Newend);
    
    /* len(): returns the length of the segment. */
    unsigned int len() const { return(End-Beg+1); }
    
    /* strand_no(): this always returns 1 since a linear segment is
     * a single contiguous stretch of the chain. Apparently meaningless
     * but nevertheless needed in virtual calls when a mixed list
     * of linear segments and sheets are processed
     */
    unsigned int strand_no() const { return(1); }
    
    /* member(): returns true if Beg<=Resno<=End, false otherwise. */
    bool member(unsigned int Resno) const { return(bool(Beg<=Resno && Resno<=End)); }

    /* mask(): constructs and returns a bitmap in which all positions
     * corresponding to the members of the segment are set true, the others
     * false. Since Linsegm_ objects know nothing about the model chain length, 
     * this should be supplied in Rno. If Rno==0 (the default) then the
     * mask will be end()+1 positions long which can be set to the correct
     * chain length later. If Rno>0 but is less than enough for containing the
     * segment then a warning is printed and an end()+1 -long mask is created.
     */
    Bits_ mask(unsigned int Rno=0) const;
};
// END OF CLASS Linsegm_

/* Class Strand_ : implements a single beta-strand. Derived from Linsegm_ .
 * Used mainly as a component of the Sheet_ class (see below).
 * A Strand_ can have a Sense (w/ respect to the previous strand in
 * the sheet a la PDB) and a Phase (also w/ respect to the previous strand).
 * Phase==0 for the first strand, Phase=Beg(this)-Beg(prev) for PAR strands, 
 * Phase=Beg(this)-End(prev) for ANTI strands.
 */
class Strand_ : public Linsegm_
{
    // typedef
    public:
    
    /* Sense_ : enum type for coding the sense of this strand w/ respect
     * to the previous strand. ANTI==-1 is antiparallel, PAR==1 is parallel, 
     * NONE==0 is for the first strand in the sheet. The int values follow
     * the PDB convention.
     */
    enum Sense_ {ANTI=-1, NONE, PAR};
    
    // data
    protected:
    
    Sense_ Sense;   // sense of strand w/ respect to previous strand
    int Phase;	    // phase of strand w/ respect to previous strand
    
    // methods
    public:
    
	// constructor
    /* Makes a Strand_ beginning at Start, ending at Stop (default==0
     * for both) with sense S (default==NONE). Phase is set to 0: will be
     * meaningful only when the object is used within a Sheet_ .
     */
    Strand_(unsigned int Start=0, unsigned int Stop=0, Sense_ S=NONE) : 
	Linsegm_(Start, Stop), Sense(S), Phase(0) {}
    
	// access
    /* sense(): w/o arguments, returns the Sense value.
     * With an int argument S, sets the sense to ANTI if S<0, 
     * to NONE if S==0, PAR if S>0 and returns the old sense value.
     */
    Sense_ sense() const { return(Sense); }
    Sense_ sense(int S);
    
    /* phase(): w/o arguments, returns the phase.
     * With an int argument P, sets the phase to P and returns the
     * old phase.
     */
    int phase() const { return(Phase); }
    int phase(int P);
};
// END OF CLASS Strand_

/* Class Sheet_ : implements a beta-sheet composed of >=2 strands.
 * Derived from Segmbase_ . Can represent intramolecular sheets only
 * and cannot handle bulges, irregularities etc.
 */
class Sheet_ : public virtual Segmbase_
{
    // data
    protected:
    
    Array_<Strand_> Strands;	// array of strands
    
    // methods
    public:
    
	// constructors
    /* Default constructor: makes a 1-strand "sheet" with the strand containing
     * the 0-th residue only. This is hardly meaningful but as sheets will be
     * built gradually, any starting point is good. Uses the default ctor of
     * Linsegm_ implicitly when the Strands array is constructed.
     */
    Sheet_() : Strands(1) {}
    
    /* Inits with a single strand. Note that 1-strand sheets do not exist
     * so use this just as a starting point for building sheets.
     * Sets the first strand's phase to 0 and sense to NONE.
     */
    Sheet_(const Strand_& Str1);
    
	// access
    /* strand_no(): returns the number of strands in the sheet. */
    unsigned int strand_no() const { return(Strands.len()); }
    
    /* [Sno]: returns the Sno-th strand of the sheet. If Sno is
     * invalid, a warning is printed and the 0-th strand is returned.
     * Note that the returned value cannot be an lvalue.
     */
    const Strand_& operator[](unsigned int Sno) const;
    
    /* member(): returns true if Resno is within the sheet, false otherwise. */
    bool member(unsigned int Resno) const;
    
    /* mask(): returns a Bits_ object in which all the bits corresponding to
     * the members of the sheet are switched ON, the rest OFF. Since Sheet_
     * objects know nothing about the size of the chain, the chain length
     * Rno must be specified. If it is omitted, then a default Rno==0
     * is assumed which means that the Bits_ result will have a length suitable
     * for accommodating all members. It can then be enlarged later.
     * If Rno>0 is specified but it is too small, then it is adjusted.
     */
    Bits_ mask(unsigned int Rno=0) const;
    
	// sheet building
    /* first_strand(): completely resets the sheet and adds Str1
     * as its first strand. The sense and phase information in Str1
     * is ignored: the phase is set to 0 and the sense to NONE inside.
     * If the sheet was initialised w/ a strand (see the relevant ctor)
     * then this method need not be used. 
     */
    void first_strand(const Strand_& Str1);
    
    /* add_strand(): adds a new strand Str to the calling object in
     * the orientation specified within, phased so that the Thisres-th
     * residue in the new strand matches the Otherres-th residue in
     * the last strand within. No action is taken if Str.sense()==NONE
     * or if Thisres and Otherres are not within their respective strands.
     * Return value: 0 on error, the new length if successful.
     */
    int add_strand(const Strand_& Str, int Thisres, int Otherres);
    
	// strand positions
    /* sense(): returns the relative sense of orientation between
     * the strands no. S1 and S2. Returns NONE if S1==S2 or if 
     * S1 or S2 are not valid strand numbers (with a warning).
     */
    Strand_::Sense_ sense(unsigned int S1, unsigned int S2) const;
	
    /* strand_res(): returns the index of the strand which contains
     * the residue Resno. Returns -1 (NOT 0!!) if Resno is not in the sheet.
     * Functionally similar to the Linsegm_::member() method.
     */
    int strand_res(unsigned int Resno) const;
    
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
    int offs_strd(unsigned int Sno, int Offs) const;
};
// END OF CLASS Sheet_

// ==== END OF HEADER Segment.h ====
#endif	/* SEGMENT_CLASSES */

