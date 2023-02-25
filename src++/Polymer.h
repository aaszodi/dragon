#ifndef POLYMER_CLASS
#define POLYMER_CLASS

// ==== PROJECT DRAGON: HEADER Polymer.h ====

/* The Polymer_ class holds all information about the model
 * chain,  viz. sequence, hydrophobicity, conservation, 
 * accessibility etc.*/

// SGI C++ 4.0, IRIX 5.3, 1. July 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "String.h"

// ---- MODULE HEADERS ----

#include "Property.h"
#include "Simil.h"
#include "Align.h"
#include "Distpred.h"
#include "Acdist.h"

// ---- IMPORTANT NOTE ----

/* The polypeptide chain is modelled as a C-alpha (CA) backbone with
 * side chains represented by one fake "C-beta" atom corresponding
 * to the side chain centroid (SCC). For traditional reasons, 
 * the SCC-s are often referred to as "fake C-betas" or simply CB-s.
 * At the termini, the model chain contains two "pseudo"-C-alpha
 * points corresponding to the NH3+ and COO- moieties. These are
 * not represented here but do play an important role in determining
 * the SCC geometry (cf. "Fakebeta" etc.). The distance matrices and
 * coordinate arrays are dimensioned so that the size is larger by 2
 * than the number of monomers returned by Polymer_::len() and 
 * the 0-th and len()+1th positions hold info about the termini.
 */

// ==== CLASSES ====

/* Class Polymer_ : keeps (almost) all information about the model chain.
 * Contains an array of Monomer_ objects for sequence, accessibility
 * etc., and internal alignment, consensus and property information.
 * Data inside a Polymer_ object can be modified via reading disk files
 * only to maintain consistency.
 */
class Polymer_
{
    /* struct Monomer_ : holds information about one amino acid in the chain.
     * Similar to the monomer struct in DRAGON-3.x but does not hold
     * secondary structure prediction information.
     */
    struct Monomer_
    {
	// data
	char Aa;	// 1-letter amino acid code
	
	double Cons;    // normalised conservation value [0..1]
	double Phob;    // average hydrophobicity
	
	double Bumpb;	    // fake C-beta bump radius (NOT squared)
	double Bumpab;	    // CA:CB bump radius (squared)
	double Abdist;   // C-alpha:sidechain centroid distance (squared)
	
	// constructor
	Monomer_(): Aa('X'), Cons(1.0), Phob(0.0), Bumpb(0.0), Bumpab(0.0), Abdist(0.0) {}
    };
    
    // static data
    static Property_ Hyphob;	// hydrophobicity
    static Property_ Volume;	// side chain volume
    static Simil_ Simil;	// amino acid similarity matrix
    static Distpred_ Dp;	// hydrophobic distance prediction
    static Acdist_ Acdist;	// sidechain atom distances from (CA|SCC)

    // data
    Align_ Align;	    // sequence alignment
    unsigned int Master;    // index of master sequence-1 (0 for consensus)
    Array_<Monomer_> Monomers;    // the polymer sequence
    Array_<double> Consphob;	// conservation*phobicity
    float Cavg, Csd;	// average and SD of phobicity
    int Changed;    // currently non-0 if parameter estimation is needed
    
    /* The polymer sequence and its properties can be changed only via
     * reading from disk files. Some changes will necessitate others and
     * the dependencies are best to deal with in a unified input routine.
     * The enum-s below tell the input routines what to modify (they may
     * be combined by | and & to achieve masking).
     */
    enum {HYPHOB=0x1, VOLUME=0x2, ACDIST=0x4, SIMIL=0x8, ALIGN=0x10};
    
    // methods
    public:
    
	// constructors
    Polymer_() : Align(), Master(0), Monomers(0), Consphob(0), 
	    Cavg(0.0), Csd(0.0), Changed(1) {}
    
	// access
    /* len(): the length of the target polymer (w/o the 0th and N+1th fake C-alphas) */
    unsigned int len() const { return(Monomers.len()); }
    
    /* cons_avg(), cons_sd(): the average and SD of the conservation
     * values in the target polymer.
     */
    float cons_avg() const { return(Cavg); }
    float cons_sd() const { return(Csd); }
    
    /* The following methods return some property of the Idx-th monomer.
     * Access is safe, indexing assumes a range of [0..Rno-1].
     * NOTE: bumpb(i) returns the NON-squared fake C-beta bump radius.
     * bumpab(i) returns the squared minimal C-alpha:C-beta bump distance.
     */
    const char& aa(unsigned int Idx) const { return(Monomers(Idx).Aa); }
    const double& cons(unsigned int Idx) const { return(Monomers(Idx).Cons); }
    const double& phob(unsigned int Idx) const { return(Monomers(Idx).Phob); }
    const double& bumpb(unsigned int Idx) const { return(Monomers(Idx).Bumpb); }
    const double& bumpab(unsigned int Idx) const { return(Monomers(Idx).Bumpab); }
    const double& abdist(unsigned int Idx) const { return(Monomers(Idx).Abdist); }

    // const access is allowed to the alignment sub-object
    const Align_& align() const { return(Align); }
    
    /* The following methods return the distance of an Atom (specified
     * as a 4-char string (PDB style)) either from the C-alpha atom
     * or from the side-chain centroid in the Idx-th residue. If Atom
     * is nonexistant in the Idx-th amino acid, then -1.0 is returned.
     * The distances are NOT squared.
     */
    float ca_dist(unsigned int Idx, const String_& Atom) const
    { return(Acdist.ca_dist(Monomers(Idx).Aa, Atom)); }

    float scc_dist(unsigned int Idx, const String_& Atom) const
    { return(Acdist.scc_dist(Monomers(Idx).Aa, Atom)); }

    /* estim_dist(): returns the NON-SQUARED estimated distance between
     * residues R1, R2 estimated from their conserved hydrophobicity
     * scores. Prints a warning and returns 0.0 if R1, R2 are out of
     * range.
     */
    double estim_dist(unsigned int R1, unsigned int R2);
    
    /* master(): returns 0 if the master sequence is the consensus of the
     * alignment, and i+1 if the i-th sequence in the alignment is the master.
     * master(Mseq): changes the master sequence within the alignment. An
     * argument of 0 means the consensus, otherwise the Mseq-1:th sequence
     * will be the consensus. No action is taken if Mseq is invalid.
     * Returns old master sequence number.
     */
    unsigned int master() const { return(Master); }
    unsigned int master(unsigned int Mseq=0);
    
    /* seq_simil(): calculates a rough sequence similarity between
     * the S1:th and S2:th sequences in the current alignment using
     * the current similarity matrix. S1, S2 in [0..Seqno]. If S1 or S2
     * are -1, then the consensus sequence of the alignment is used.
     * The score returned is an average of pairwise scores between all non-gapped
     * positions. Returns -1.0 on error.
     */
    float seq_simil(int S1, int S2) const;
    
	// input
    /* NOTE: A Polymer_ object has several components with individual
     * input routines. Overloading of >> cannot be done. The read_*()
     * methods read from files, the str_*() methods from strings.
     */
    
    /* read_aln(): attempts to read an alignment file from Fname. If Mseq==0,
     * then the "master sequence" will be the consensus of the alignment (default);
     * otherwise, the (Mseq-1)-th sequence from the alignment will be the
     * master sequence (i.e. the sequence of the model chain). If the input
     * is not successful or there are not enough sequences then nothing will
     * be changed; however, a successful alignment input will cause all the
     * other fields be updated.
     * Return value: 0 on failure, the new sequence length on success.
     */
    int read_aln(const char *Fname, unsigned int Mseq=0);
    
    /* str_aln(): does the same thing as read_aln() but reads from 
     * a string Str.
     */
    int str_aln(const char *Str, unsigned int Mseq=0);
    
    /* read_phob(), read_vol(), read_acdist, read_simil(): 
     * read new data from a file Fname. Update the
     * sequence records accordingly. Return 0 on failure, >=1 on success.
     * The corresponding str_*() methods do exactly the same but
     * read the input from a string Str rather than from a file.
     */
    int read_phob(const char *Fname);
    int str_phob(const char *Str);
    int read_vol(const char *Fname);
    int str_vol(const char *Str);
    int read_acdist(const char *Fname);
    int str_acdist(const char *Str);
    int read_simil(const char *Fname);
    int str_simil(const char *Str);
    
	// output
    friend ostream& operator<<(ostream& Out, const Polymer_& P);
    
	// private input
    private:
    void update_members(int Mod);
    
	// "forbidden methods": no copy/assignment
    Polymer_(const Polymer_&);
    Polymer_& operator=(const Polymer_&);
};
// END OF CLASS Polymer_

// ==== END OF HEADER Polymer.h ====

#endif	/* POLYMER_CLASS */
