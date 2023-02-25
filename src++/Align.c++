// ==== PROJECT DRAGON: METHODS Align.c++ ====

/* Class for storing multiple alignments. */

// SGI C++, 12-Feb-1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <strstream.h>
#include <ctype.h>
#include <string.h>

// ---- CLASS HEADER ----

#include "Align.h"

// ---- STATIC DEFINITIONS ----

const int Align_::MAXSEQLEN=2048;
const int Align_::MAXSEQNO=256;
const char Align_::GAP='-';
const int Align_::LINELEN=81;

// ==== Align_ METHODS ====

// ---- Constructors ----

/* Inits to hold MAXSEQLEN strings by default: members are NULL. */
Align_::Align_(): Len(MAXSEQLEN), Seqno(0)
{
    Aln=new  char* [MAXSEQLEN];
    for (register unsigned int i=0; i<MAXSEQLEN; i++) Aln[i]=NULL;
}

// ---- Destructor ----

Align_::~Align_()
{
    for (register unsigned int i=0; i<Len; i++) delete [] Aln[i];
    delete [] Aln;
}

// ---- Access ----

/* seq(): returns the Idx-th sequence in Seq the size of which will be
 * appropriately adjusted. If Idx is out of range then a warning is
 * printed, nothing is done and 0 is returned,  otherwise the 
 * actual sequence length is returned. This can be shorter 
 * than the safely allocated array because the gaps are skipped.
 */
unsigned int Align_::seq(unsigned int Idx, char *& Seq) const
{
    if (Idx>=Seqno)
    {
	cerr<<"\n? Align_::seq("<<Idx<<"): Out of range\n";
	return(0);
    }
    
    if (Seq==NULL || strlen(Seq)<=Len)	// realloc
    {
	delete [] Seq;	    // does nothing with NULL
	Seq=new char [Len+1];
    }
    memset(Seq, '\0', Len+1);	// paranoia
    
    register unsigned int i, p;
    register char Ct;
    for (i=p=0; p<Len; p++)	// extract Idx-th col
    {
	Ct=Aln[p][Idx];
	if (Ct==GAP) continue;
	else Seq[i++]=Ct;
    }
	
    return(i);    // OK
}
// END of seq()

/* seq_len(): returns the "net length" of the Idx-th sequence in the
 * alignment, i.e. the number of positions minus the gaps.
 */
unsigned int Align_::seq_len(unsigned int Idx) const
{
    if (Idx>=Seqno)
    {
	cerr<<"\n? Align_::seq_len("<<Idx<<"): Out of range\n";
	return(0);
    }
    
    register unsigned int i, L=Len; // overall length
    for (i=0; i<Len; i++)
	if (Aln[i][Idx]==GAP) L--;
    return(L);
}
// END of seq_len()

/* pos(): returns the Idx-th alignment position as a const char*.
 * If Idx is out of range then a warning is printed and NULL is returned.
 */
const char* Align_::pos(unsigned int Idx) const
{
    if (!Seqno)
    {
	cerr<<"\n? Align_::pos(): No sequences\n";
	return(NULL);
    }
    
    if (Idx>=Len)
    {
	cerr<<"\n? Align_::pos("<<Idx<<"): Out of range\n";
	return(NULL);
    }
    
    return(Aln[Idx]);  // OK
}
// END of pos()

/* align_pos(): returns the alignment position which contains
 * the Pos:th position of the Idx:th sequence. Idx must be 
 * in the range [0..Seqno-1]. Pos is smaller than the length of 
 * the sequence. The returned position falls within [0..Len-1], 
 * -1 or -2 indicates that Idx or Pos was invalid.
 */
int Align_::align_pos(unsigned int Idx, unsigned int Pos) const
{
    if (Idx>=Seqno)
    {
	cerr<<"\n? Align_::align_pos("<<Idx<<", ...): Out of range\n";
	return(-1);
    }
    if (Pos>=Len)
    {
	cerr<<"\n? Align_::align_pos(..., "<<Pos<<"): Out of range\n";
	return(-2);
    }
    
    register unsigned int i, p;
    for (i=p=0; p<=Pos && i<Len; i++)
    {
	if (Aln[i][Idx]==GAP) continue;
	++p;
    }
    return(i-1);
}
// END of align_pos()

/* seq_pos(): given the alignment position Pos, the corresponding
 * sequence position of the Idx:th sequence is returned. If there
 * is a gap, then -1 is returned. Invalid ranges result in warnings
 * and -2 is returned.
 */
int Align_::seq_pos(unsigned int Idx, unsigned int Pos) const
{
    if (Idx>=Seqno)
    {
	cerr<<"\n? Align_::seq_pos("<<Idx<<", ...): Out of range\n";
	return(-2);
    }
    if (Pos>=Len)
    {
	cerr<<"\n? Align_::seq_pos(..., "<<Pos<<"): Out of range\n";
	return(-2);
    }
    if (Aln[Pos][Idx]==GAP)	// 19-Jan-96
	return(-1);
    
    register unsigned int p, Gapno;
    for (p=Gapno=0; p<=Pos; p++)
	if (Aln[p][Idx]==GAP) ++Gapno;
    return(Pos-Gapno);
}
// END of seq_pos()

/* reset(): clears the calling object to its nascent state.
 * The Aln array will be L (default MAXSEQLEN) long and all items will be NULL.
 * If L==0 or >MAXSEQLEN then it'll be reset silently to MAXSEQLEN.
 * Seqno will always be set to 0.
 */
void Align_::reset(unsigned int L)
{
    register unsigned int i;
    for (i=0; i<Len; i++) delete [] Aln[i]; // equivalent to ~Align_()
    delete [] Aln;
    
    if (!L || L>=MAXSEQLEN) L=MAXSEQLEN;    // reset size silently if L is silly
    Aln=new char* [Len=L];
    for (i=0; i<Len; i++) Aln[i]=NULL;
    Seqno=0;
}
// END of reset()

// ---- Input ----

/* read_file(): reads a multiple alignment file from Fname.
 * Reading continues up to
 * EOF or until MAXSEQLEN positions have been read.
 * The calling object will be modified only if the whole operation
 * was successful in which case the overall length is returned: otherwise, 
 * 0 is returned.
 */
unsigned int Align_::read_file(const char *Fname)
{
    ifstream Infile(Fname);
    if (!Infile)
    {
	cerr<<"\n? Align_::read_file("<<Fname<<"): Cannot open\n";
	return(0);
    }
    Infile>>(*this);	// input in one go
    Infile.close();
    
    return((Infile.good() || Infile.eof())? len(): 0);
}
// END of read_file()

/* >>: tries to input an alignment file from Inf into Align. 
 * Tries the vertical MULTAL format, GCG's MSF, PIR format.
 * If none succeeds, then Inf's failbit is set, the stream is rewound to
 * its position prior to the input and Align is not modified.
 */
istream& operator>>(istream& Inf, Align_& Align)
{
    cout<<"\n# >> Align_: Trying MULTAL format...\n";
    if (Align.read_multal(Inf))
    {
	cout<<"# >> Align_: MULTAL parsing successful, seqno="<<Align.seq_no()<<endl;
	return(Inf);
    }
    else cout<<"# >> Align_: not in MULTAL format...\n";
    
    // not MULTAL, read_multal() has rewound Inf
    cout<<"\n# >> Align_: Trying GCG-MSF format...\n";
    if (Align.read_msf(Inf))
    {
	cout<<"# >> Align_: MSF parsing successful, seqno="<<Align.seq_no()<<endl;
	return(Inf);
    }
    else cout<<"# >> Align_: not in MSF format...\n";
    
    // try PIR format
    cout<<"\n# >> Align_: Trying PIR format...\n";
    if (Align.read_pir(Inf))
    {
	cout<<"# >> Align_: PIR parsing successful, seqno="<<Align.seq_no()<<endl;
	return(Inf);
    }
    else cout<<"# >> Align_: not in PIR format...\n";

    // bad luck
    cerr<<"\n? >>Align_: Sorry, cannot parse alignment file\n";
    Inf.clear(Inf.rdstate()|ios::failbit);
    return(Inf);
}
// END of >>

/* read_multal(): reads a MULTAL-like multiple alignment file from the stream Inf.
 * If the input was successful, then the number of sequences (>0)
 * is returned. If the file could not be parsed, then 0 is returned
 * and Inf is reset to its status before the call. Private
 */
int Align_::read_multal(istream& Inf)
{
    if (!Inf)
    {
	cerr<<"\n? Align_::read_multal(): Cannot read input stream\n";
	return(0);
    }
    streampos Origpos=Inf.tellg();   // save original position
    long Origflags=Inf.flags();	    // and format state
    int Origstate=Inf.rdstate();    // and error state

    char Line[MAXSEQNO+3];  // input buffer
    istrstream Inline(Line, MAXSEQNO+2); // space for trailing \n\0
    char *Seqptr=NULL;	// ptr to sequence name
    unsigned int Lineno, Newlen=0;
    int Sno=0, Blockseen=0;
    char **Temp=NULL;
    
    // process line by line
    for (Lineno=0; Inf && Newlen<MAXSEQLEN; Lineno++)
    {
	Line[0]='\0';	// reset input buffer line
	Inf.getline(Line, MAXSEQNO+2, '\n');
	if (Line[0]=='#' || Line[0]=='\0' || Line[0]=='\n')
	    continue;	// skip comments and empty lines
	
	// if a sequence name line is found, echo to stdout
	Seqptr=strstr(Line, "USER>");
	if (Seqptr!=NULL)
	{
	    cout<<"# Sequence:"<<(Seqptr+5)<<endl;
	    continue;
	}
	
	/* Get the number of sequences. This is somewhat difficult
	 * because there are a number of MULTAL formats around :-(
	 * DRAGON: expects "Seqno %d" at the beginning of the line
	 * CAMELEON's MULTAL format: "[bB]lock %d = %d seqs"
	 * MSAP format: "%d seqs" (preceded by a [bB]lock %d\n")
	 */
	if (!Sno)   // hasn't seen number of sequences yet
	{
	    if (strlen(Line)>=7 && !strncmp(Line, "Seqno ", 6)) // DRAGON
		Inline.seekg(6, ios::beg);  // position after "Seqno"
	    else 
	    {
		// check for the presence of the block number
		if (!Blockseen)
		    Blockseen=(Line[0]=='B' || Line[0]=='b') &&
			!strncmp(Line+1, "lock", 4);
		if (!Blockseen) continue;
		
		// assuming CAMELEON variant, find '='
		Seqptr=strchr(Line+6, '=');
		if (Seqptr==NULL) Seqptr=Line-1;    // MSAP variant?
		if (NULL!=strstr(Line, "seqs"))
		    Inline.seekg(Seqptr-Line+1, ios::beg);	// pos before "%d seqs"
		else continue;
	    }
	    
	    Inline>>Sno;
	    if (!Inline || !Sno)
	    {
		cerr<<"\n? Align_::read_multal(): Sequence no. "<<Sno<<" invalid or missing: stopped reading\n";
		Sno=0; break;
	    }
	    Sno=abs(Sno);
	    
	    // set up temp array
	    Temp=new char* [MAXSEQLEN];
	    continue;
	}	// if (!Sno)
	
	// attempt to put into temp
	if (NULL==(Temp[Newlen]=check_vertical(Line, Sno)))
	{
	    cerr<<"\n? Align_::read_multal(): Line "<<(Lineno+1)<<" cannot be parsed\n";
	    continue;
	}
	Newlen++;
    }	    // for
    
    if (!Sno)	// probably not MULTAL format
    {
	cerr<<"\n? Align_::read_multal(): Input file is not in MULTAL format\n";
	Inf.flags(Origflags); Inf.seekg(Origpos);	// reset Inf to status before call
	Inf.clear(Origstate);
	return(0);
    }
    if (Sno && Newlen)	// modify calling object
    {
	reset(Newlen);
	for (unsigned int i=0; i<Newlen; i++)
	    Aln[i]=Temp[i];
	Seqno=Sno;
    }
    delete [] Temp;
    return(Sno);    // OK
}
// END of read_multal()

/* read_msf(): attempts to read the input stream Inf into the calling object,
 * assuming it contains a multiple alignment in GCG's MSF format.
 * Returns the no. of sequences if the input operation was successful, 
 * returns 0 if Inf does not correspond to MSF format. In this case, 
 * Inf is reset to its original status before the call. Private
 */
int Align_::read_msf(istream& Inf)
{
    if (!Inf)
    {
	cerr<<"\n? Align_::read_msf(): Cannot read input stream\n";
	return(0);
    }
    streampos Origpos=Inf.tellg();   // save original position
    long Origflags=Inf.flags();	    // and format state
    int Origstate=Inf.rdstate();    // and error status
    
    char Line[MAXSEQLEN+LINELEN], Namestr[6], Names[MAXSEQNO][LINELEN], 
	Namebuf[LINELEN], Lenstr[5];  // input buffers
    istrstream Inline(Line, MAXSEQLEN+LINELEN);
    unsigned int Lineno, p, Seqlen, Namelen, Maxlen=0, Maxnamelen=0, 
	Sno=0, Newlen=0, Scur=0, Pcur=0, Err=1;
    char **Temp=NULL;
    
    // process line by line
    for (Lineno=0; Inf && Pcur<=Newlen; Lineno++)
    {
	Line[0]='\0';	// reset input buffer line
	Inf.getline(Line, MAXSEQLEN+LINELEN, '\n');
	if (Line[0]=='#' || Line[0]=='\0' || Line[0]=='\n')
	    continue;	// skip comments and empty lines
	Inline.seekg(0); Inline.clear(0);   // reset
	
	// hunt for sequence names: "Name: %s Len: %d", store name and seq lengths
	if (!Newlen)
	{
	    Inline>>setw(6)>>Namestr>>
		setw(LINELEN)>>Names[Sno]>>setw(5)>>Lenstr>>Seqlen;
	    if (!Inline || strcmp(Namestr, "Name:") || strcmp(Lenstr, "Len:"))
	    {
		// the current line is not a name line
		if (Maxlen)  // finish name processing
		{
		    Newlen=Maxlen;
		    Temp=new char* [Newlen];	// set up temporary storage
		    for (p=0; p<Newlen; p++)
		    {
			Temp[p]=new char [Sno+1];	// one position
			Temp[p][Sno]='\0';
			memset(Temp[p], GAP, Sno);  // set to "---...-"
		    }
		    Scur=Pcur=Maxlen=0; // set current positions and seq length to 0
		    Err=0;  // clear error status
		    Inline.clear(0); Inline.seekg(0);	// reset for re-processing line
		}
		else continue;	// no name line seen yet
	    }
	    else    // name line found
	    {
		cout<<"# Sequence:"<<Names[Sno]<<endl;   // echo name to stderr
		if (Seqlen>Maxlen) Maxlen=Seqlen;   // update alignment length
		Namelen=strlen(Names[Sno]);
		if (Namelen>Maxnamelen) Maxnamelen=Namelen; // sequence name length
		Sno++; continue;    // count sequences
	    }
	}
	
	/* Alignment lines are prepended with the (right-justified)
	 * sequence names. Scur indicates the current sequence, Pcur
	 * the alignment position. We allow any garbage between alignment
	 * blocks (i.e. when Scur==0), but not between lines.
	 */
	Inline>>setw(LINELEN)>>Namebuf;
	if (!Inline || !Scur && strcmp(Names[0], Namebuf)) 
	    continue;	// was not an alignment line, skip
	if (strcmp(Names[Scur], Namebuf))   // mismatch!
	{
	    cerr<<"\n? Align_::read_msf(): Seqname tag \""<<Names[Scur]
		<<"\" expected, \""<<Namebuf<<"\" found in line "<<(Lineno+1)<<endl;
	    Err=1; break;   // this is serious!
	}
	
	// tag was OK, process everything after it as alignment
	Seqlen=check_horizontal(Line+Maxnamelen);
	if (Seqlen>Maxlen) Maxlen=Seqlen;   // might have different lengths?
	for (p=Pcur; p<Pcur+Seqlen && p<Newlen; p++)
	    Temp[p][Scur]=Line[Maxnamelen+p-Pcur];
	Scur=(Scur+1)%Sno;  // increment Scur with wraparound at Sno
	if (!Scur)
	{
	    Pcur+=Maxlen;    // next chunk will come
	    Maxlen=0;	// collect chunk length here
	}
    }	    // for
    
    if (Pcur<Newlen)	// too short (this is not an error condition though)
    {
	cerr<<"\n? Align_::read_msf(): Actual alignment length is "<<Pcur
	    <<", expected "<<Newlen<<endl;
	for (p=Pcur; p<Newlen; p++) delete [] Temp[p];
	Newlen=Pcur;
    }
    if (Err || !Sno)	//  error
    {
	cerr<<"\n? Align_::read_msf(): Input file is not in MSF format\n";
	if (Temp!=NULL)
	{
	    for (p=0; p<Maxlen; p++) delete [] Temp[p];
	    delete [] Temp;
	}
	Inf.flags(Origflags); Inf.seekg(Origpos);	// reset Inf to status before call
	Inf.clear(Origstate);
	return(0);
    }

    // OK
    reset(Newlen);
    for (p=0; p<Newlen; p++) Aln[p]=Temp[p];
    Seqno=Sno;
    delete [] Temp;
    return(Sno);    // OK
}
// END of read_msf()

/* read_pir(): attempts to read the input stream Inf into the calling object,
 * assuming it contains a multiple alignment in PIR format, 
 * ie. sequences with optional gaps as PIR entries following each other.
 * This is the alignment format MODELER likes.
 * Returns the no. of sequences if the input operation was successful, 
 * returns 0 if Inf does not correspond to the PIR format. In this case, 
 * Inf is reset to its original status before the call. Private
 */
int Align_::read_pir(istream& Inf)
{
    if (!Inf)
    {
	cerr<<"\n? Align_::read_pir(): Cannot read input stream\n";
	return(0);
    }
    streampos Origpos=Inf.tellg();   // save original position
    long Origflags=Inf.flags();	    // and format state
    int Origstate=Inf.rdstate();    // and error status
    
    char Line[MAXSEQLEN+LINELEN];
    char **Temp=NULL;    // temporary storage
    char *Cptr;
    int Lineno, Sno=0, s, Chunklen, Pcur=0, p, Maxlen=0;
    
    /* every PIR entry has a "P1 line" which looks like
     * ">P1;seq_name"
     * then a second comment line which is ignored, 
     * then the sequence lines, terminated by a "*".
     * Pirstatus stores what to parse next.
     */
    enum { P1_LINE, SECOND_LINE, SEQ_LINE } Pirstatus=P1_LINE;
    
    // make room for MAXSEQNO sequences
    Temp=new char* [MAXSEQNO];
    
    // process line by line
    for (Lineno=0; Inf; Lineno++)
    {
	Line[0]='\0';	// reset input buffer line
	Inf.getline(Line, MAXSEQLEN+LINELEN, '\n');
	if (Line[0]=='#' || Line[0]=='\0' || Line[0]=='\n')
	    continue;	// skip comments and empty lines
	
	// depending on Pirstatus, decide what to parse
	switch (Pirstatus)
	{
	    case P1_LINE:   // look for a new sequence name
	    Cptr=strstr(Line, ">P1;");
	    if (Cptr!=NULL)	// first line found
	    {
		cout<<"# Sequence: "<<(Cptr+4)<<endl;	// print name
		Pirstatus=SECOND_LINE;	// now look for second line
	    }
	    break;
	    
	    case SECOND_LINE:	// this is just a comment
	    cout<<"# Description: "<<Line<<endl;
	    Pirstatus=SEQ_LINE;	    // now look for the (Sno:th) sequence itself
	    Temp[Sno]= new char [MAXSEQLEN+1];	// make room for it
	    Pcur=0; // init current position in temp storage
	    break;
	    
	    case SEQ_LINE:	// process sequence info
	    Cptr=strrchr(Line, '*');	// locate terminating asterisk
	    if (Cptr!=NULL)
	    {
		*Cptr='\0';	// terminate chunk at asterisk
		Pirstatus=P1_LINE;  // and look for new sequence
	    }
	    Chunklen=check_horizontal(Line);
	    if (Pcur+Chunklen>MAXSEQLEN)
	    {
		cerr<<"\n? >>Align_::read_pir(): line "<<
		    (Lineno+1)<<": Sequence too long ('*' missing?)\n";
		Pirstatus=P1_LINE;
		break;	// throw whole line away: perhaps too drastic
	    }
	    strcpy(Temp[Sno]+Pcur, Line);   // transfer cleansed chunk
	    Pcur+=Chunklen;
	    if (Pirstatus==P1_LINE)
	    {
		if (Pcur>Maxlen) Maxlen=Pcur;	// save maximal sequence length
		Sno++;  // prepare for next round
	    }
	    break;
	}	// switch Pirstatus
    }	    // for Lineno
    
    if (!Sno || !Maxlen)    // something went wrong
    {
	cerr<<"\n? Align_::read_pir(): Input file is not in PIR format\n";
	Inf.flags(Origflags); Inf.seekg(Origpos);	// reset Inf to status before call
	Inf.clear(Origstate);
	Sno=0;
    }
    else	// was OK
    {
	// check if all sequences have the same length, pad shorter ones with gaps
	int Slen;
	for (s=0; s<Sno; s++)
	{
	    Slen=strlen(Temp[s]);
	    if (Slen<Maxlen)
	    {
		memset(Temp[s]+Slen, '-', Maxlen-Slen);
		Temp[s][Maxlen]='\0';
		cerr<<"\n? Align_::read_pir(): Sequence "<<(s+1)
		    <<" too short ("<<Slen<<"<"<<Maxlen<<"), padded with gaps\n";
	    }
	}
	// transfer into vertical storage
	reset(Maxlen);
	for (p=0; p<Maxlen; p++)
	{
	    Aln[p]=new char [Sno+1];
	    for (s=0; s<Sno; s++) Aln[p][s]=Temp[s][p];
	    Aln[p][Sno]='\0';
	}
	Seqno=Sno;
    }
    
    // cleanup
    if (Temp!=NULL)
    {
	for (s=0; s<Sno; s++) delete [] Temp[s];
	delete [] Temp;
    }
    return(Sno);
}
// END of read_pir() 

/* check_vertical(): checks if the input string Instr corresponds to
 * the vertical MULTAL alignment format. Should be Sno long, must contain at least one
 * alpha character. Lowercase alphas are converted to uppercase, anything
 * else will be converted to GAP signs '-'. Allocates and returns a
 * Sno-long string (padded at the end with '-'s if necessary) or
 * returns NULL if Instr was really hopeless (i.e. zero-long or 
 * gap signs only). Private static
 */
char* Align_::check_vertical(char *Instr, int Sno)
{
    unsigned int i, Inlen, Ano=0;
    char *Nlpos=NULL;
    char Aa;
    
    Nlpos=strrchr(Instr, '\n');	// locate terminal \n
    if (Nlpos!=NULL) *Nlpos='\0';
    Inlen=strlen(Instr);    // now get the length
    
    for (i=0; i<Inlen; i++)
    {
	Aa=Instr[i];
	if (isalpha(Aa)) { Instr[i]=toupper(Aa); Ano++; }
	else if (Aa!=GAP)
	{
	    cerr<<"\n? Align_::check_vertical("<<Instr<<"): Illegal AA code '"<<Aa<<"', replaced by 'X'\n";
	    Instr[i]='X'; Ano++;
	}
    }
    if (!Ano) return(NULL); // hopeless: not an alpha!
    
    // compare the lengths
    if (Inlen!=Sno)
    {
	cerr<<"\n? Align_::check_vertical("<<Instr<<"): No. of positions "<<Inlen;
	if (Inlen<Sno)
	    cerr<<"<"<<Sno<<", padded with '-'\n";
	else
	    cerr<<">"<<Sno<<", truncated\n";
    }
    
    char *Outstr=new char [Sno+1];
    Outstr[Sno]='\0'; memset(Outstr, GAP, Sno);	    // fill up with all '-'
    strncpy(Outstr, Instr, (Inlen>Sno)? Sno: Inlen);	// trunc if necessary
    return(Outstr);
}
// END of check_vertical()

/* check_horizontal(): checks if the argument string Instr corresponds to
 * the horizontal alignment format. Currently this is GCG's MSF, but
 * the usual gap character '-' is also accepted. Instr is modified
 * as follows: lowercase characters are converted to uppercase, the MSF
 * gap characters '.' are changed to '-', unrecognised characters are
 * replaced by 'X' with a warning, whitespaces removed.
 * Return value: the new length of Instr. Private static
 */
int Align_::check_horizontal(char *Instr)
{
    int p, q, Inlen=strlen(Instr);
    char Aa;
    
    for (p=q=0; p<Inlen; p++)
    {
	Aa=Instr[p];
	if (isalpha(Aa)) { Instr[q++]=toupper(Aa); continue; }	// [a-z]-->[A-Z]
	switch(Aa)
	{
	    case GAP: case '.': Instr[q++]=GAP; break;	// MSF or normal gap-> normal gap
	    case ' ': case '\t': case '\n': break;  // skip ws
	    default:
		cerr<<"\n? Align_::check_horizontal(): Illegal AA code \'"<<
		    Aa<<"\', replaced by \'X\'\n";
		Instr[q++]='X';
	    break;
	}
    }
    Instr[q]='\0';  // terminate and return
    return(q);
}
// END of check_horizontal()

/* <<: output to ostream Out, in a very simple format,
 * just the sequences under each other, no extra info.
 */
ostream& operator<<(ostream& Out, const Align_& A)
{
    static const int SPACE_INTERVAL=10, CHARS_PER_LINE=60;
    int p, pst, s;
    
    for (pst=0; pst<A.len(); pst+=CHARS_PER_LINE)
    {
	for (s=0; s<A.seq_no(); s++)
	{
	    for (p=pst; p<pst+CHARS_PER_LINE && p<A.len(); p++)
	    {
		Out<<A.Aln[p][s];
		if ((p+1)%SPACE_INTERVAL ==0) Out<<' ';
	    }
	    Out<<endl;
	}
	Out<<endl;
    }
    return(Out);
}
// END of <<

// ==== END OF METHODS Align.c++ ====

