#ifndef ALIGN_CLASS
#define ALIGN_CLASS

// ==== PROJECT DRAGON: HEADER Align.h ====

/* Class for storing multiple alignments. */

// SGI C++, 12-Feb-1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>

// ==== CLASSES ====

/* Class Align_: stores a multiple alignment. Modifiable only
 * by the read_file() method that understands a MULTAL-like
 * input. Can be queried for individual sequences and
 * individual alignment positions.
 */
class Align_
{
    // data
    private:

    static const int MAXSEQLEN;	// maximal sequence length
    static const int MAXSEQNO;	// maximal number of sequences
    static const char GAP;  // the gap character
    static const int LINELEN;	// the maximal width of a normal input line

    char **Aln;	// like in MULTAL, sequences run "vertically" (Len x Seqno array)
    unsigned int Len, Seqno;	// Seqno sequences, total length is Len
    
    // methods
    public:
    
	// constructors
    /* Inits to hold MAXSEQLEN strings by default: members are NULL. */
    Align_();
    
	// destructor
    ~Align_();
    
	// access
    
    /* The sizes can be changed by read_file() only. Note that
     * len() returns 0 if there are no sequences even though the
     * allocated length of the alignment is MAXSEQLEN inside
     */
    unsigned int len() const { return(Seqno? Len: 0); }
    unsigned int seq_no() const { return(Seqno); }
    
    /* seq(): returns the Idx-th sequence in Seq the size of which will be
     * appropriately adjusted. If Idx is out of range then a warning is
     * printed, nothing is done and 0 is returned,  otherwise the 
     * actual sequence length is returned. This can be shorter 
     * than the safely allocated array because the gaps are skipped.
     */
    unsigned int seq(unsigned int Idx, char *& Seq) const;
    
    /* seq_len(): returns the "net length" of the Idx-th sequence in the
     * alignment, i.e. the number of positions minus the gaps.
     */
    unsigned int seq_len(unsigned int Idx) const;
    
    /* pos(): returns the Idx-th alignment position as a const char*.
     * If Idx is out of range then a warning is printed and NULL is returned.
     */
    const char* pos(unsigned int Idx) const;
    
    /* align_pos(): returns the alignment position which contains
     * the Pos:th position of the Idx:th sequence. Idx must be 
     * in the range [0..Seqno-1]. Pos is smaller than the length of 
     * the sequence. The returned position falls within [0..Len-1], 
     * -1 or -2 indicates that Idx or Pos was invalid.
     */
    int align_pos(unsigned int Idx, unsigned int Pos) const;
    
    /* seq_pos(): given the alignment position Pos, the corresponding
     * sequence position of the Idx:th sequence is returned. If there
     * is a gap, then -1 is returned. Invalid ranges result in warnings
     * and -2 is returned.
     */
    int seq_pos(unsigned int Idx, unsigned int Pos) const;
    
    /* reset(): clears the calling object to its nascent state.
     * The Aln array will be L (default MAXSEQLEN) long and all items will be NULL.
     * If L==0 or >MAXSEQLEN then it'll be reset silently to MAXSEQLEN.
     * Seqno will always be set to 0.
     */
    void reset(unsigned int L=MAXSEQLEN);
    
	// input
    
    /* NOTE: The following formats are supported:-
     * 
     * 1) MULTAL-like vertical format:
     * Perfect MULTAL compatibility cannot be achieved, due to the various
     * output formats in circulation. The format recognised here is fairly
     * flexible. The basic trick is to get the number of sequences from
     * the file. The first non-empty, non-comment line(s) of the file 
     * (the header) should look like as one of the following possibilities
     * (using scanf-like format syntax for simplicity):
     * 
     * "Seqno %d"  - the DRAGON format
     * or 
     * "Block 0\n%d seqs" - the MSAP format (capital 'B' not significant)
     * or
     * "block 1 = %d seqs" - the CAMELEON/MULTAL format
     * 
     * where %d is the number of sequences (negative numbers will be
     * converted to abs values, 0 is an error). In all cases no leading
     * whitespaces are allowed in the header.
     * The names of the sequences may be specified on the following
     * lines after a "USER>" token.
     * The sequence lines each should contain XX characters from the set
     * [a-zA-Z] or '-' for gaps. No whitespaces are allowed. If there are
     * less than XX characters, the rest is padded with '-': lines longer
     * than XX chars are truncated. Lines containing
     * only gap signs "---" are skipped w/ warnings. Illegal characters
     * trigger a warning and will be replaced by 'X'.
     * 
     * 2) MSF (Multiple Sequence Format)-like horizontal format:
     * This format is used by the GCG suite of programs. The parser
     * requires the presence of the "Name" lines:
     * "Name: %s Len: %d"
     * where the %s is the name of the sequence and %d is the length.
     * The number of sequences will be deduced from the number of "Name"
     * lines, the total length of the alignment is assumed to be the maximal
     * Len value (if they're not equal). The alignment lines should
     * look like:
     * 
     * name1  Al....IG n.ME...NT...
     * name2  aL....VG N-ME---Nt...
     * ...
     * 
     * i.e. lower-and uppercases are accepted, spaces are ignored, 
     * and both '.' (the MSF gap) and '-' (the MULTAL/DRAGON gap)
     * can be used. The sequence names (right-justified) must correspond to the
     * ones given in the Name lines in the same order.
     * 
     * 3) PIR horizontal format (Modeler input):
     * The alignment file contains one or more PIR sequences with gaps.
     * A PIR sequence looks like this:
     * 
     * >P1;sequence_name
     * description_line
     * first_sequence_line...
     * ...
     * ...last_sequence_line*
     * 
     * terminated with an asterisk.
     */
    
    /* read_file(): reads a multiple alignment file from Fname.
     * Reading continues up to
     * EOF or until MAXSEQLEN positions have been read.
     * The calling object will be modified only if the whole operation
     * was successful in which case the overall length is returned: otherwise, 
     * 0 is returned.
     */
    unsigned int read_file(const char *Fname);
    
    /* >>: tries to input an alignment file from Inf into Align. 
     * Tries the vertical MULTAL format, GCG's MSF, PIR format.
     * If none succeeds, then Inf's failbit is set, the stream is rewound to
     * its position prior to the input and Align is not modified.
     */
    friend istream& operator>>(istream& Inf, Align_& Align);
    
    /* <<: output to ostream Out, in a very simple format,
     * just the sequences under each other, no extra info.
     */
    friend ostream& operator<<(ostream& Out, const Align_& A);
    
    private:
    int read_multal(istream& Inf);
    int read_msf(istream& Inf);
    int read_pir(istream& Inf);
    static char* check_vertical(char *Instr, int Sno);
    static int check_horizontal(char *Instr);

    // "forbidden methods": no copying or assignment
    Align_(const Align_&);
    Align_& operator=(const Align_&);
};
// END OF CLASS Align_

// ==== END OF HEADER Align.h ====
#endif	/* ALIGN_CLASS */
