// ==== PROJECT DRAGON: METHODS Simil.c++ ====

/* Implements a class for amino acid similarity matrices. */

// SGI C++, IRIX 6.2, 22. Nov. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <fstream.h>
#include <strstream.h>

// ---- CLASS HEADER ----

#include "Simil.h"

// ---- DEFINITIONS ----

#define LINELEN 240

// ==== Simil_ METHODS ====

// ---- Constructors ----

/* Inits to the 26 uppercase letters of the alphabet and to the unit matrix. */
Simil_::Simil_(): Sim(26)
{
    Aacodes=new char [27];
    strcpy(Aacodes, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    Sim.set_values(); Sim.diag_matrix();    // unit matrix
}

// ---- Access ----

/* simil(): returns the similarity value between amino acids A and B.
 * Lowercase letters will be silently converted to uppercase, unknown
 * characters trigger a warning and will be converted to 'X'. If X is
 * not in the code string then 0.0 is returned.
 */
double Simil_::simil(char A, char B) const
{
    int Ai=pos(A), Bi=pos(B);
    if (Ai>=0 && Bi>=0) return(Sim(Ai, Bi)); else return(0.0);
}
// END of simil();

/* cons(): calculates the normalised consensus value for the amino acid
 * codes in the string Aas. Empty strings have a value of 0.0, strings
 * with just 1 char (possibly coming from 1-sequence alignments) give
 * 1.0, "AAA..A" have also 1.0, anybody else between 0 and 1. Gaps are
 * skipped in the summation (indicated by '-'). Returns the consensus
 * character: the most abundant amino acid in the string. The consensus
 * value is returned in Consval. (The consensus char is not always needed.)
 */
char Simil_::cons(const char *Aas, double& Consval) const
{
    // no amino acids specified
    if (Aas==NULL || Aas[0]=='\0')
    {
	Consval=0.0; return('X');
    }
    
    unsigned int Aano=strlen(Aas);
    const char GAP='-';

    // only one character in string: total consensus (even if nonsense!)
    if (Aano==1)
    {
	Consval=(Aas[0]!=GAP)? 1.0: 0.0;
        return(Aas[0]);
    }
    
    register unsigned int j, k, Pno=(Aano*(Aano-1))/2;	// pair no
    register int jx, kx;
    double Sco, Maxsco=-9.9e10;
    const unsigned int Aacodelen=strlen(Aacodes);   // code string length
    
    // prepare for counting characters in the string
    int *Anum=new int [Aacodelen];    // one pos for each char in the code string
    memset(Anum, 0, Aacodelen*sizeof(int));   // zero it
    
    Consval=0.0;
    char Aj, Ak;
    for (j=0; j<Aano; j++)  // outer cycle: visits all positions once
    {
	Aj=Aas[j];
	if (Aj==GAP || 0>(jx=pos(Aj))) continue;    // skip nonsense
	
	// get maximal possible score
	Sco=simil(Aj, Aj);
	if (Sco>Maxsco) Maxsco=Sco;

	// collect pairwise similarity scores
	for (k=0; k<j; k++)
	{
	    Ak=Aas[k];
	    if (Ak==GAP || 0>(kx=pos(Ak))) continue;
	    Sco=simil(Aj, Ak);
	    Consval+=Sco;
	    if (Sco>Maxsco) Maxsco=Sco;	// this is paranoid: self-sim is always higher than inter-sim
	}
	
	Anum[jx]++; // count the amino acids
    }
    
    // get normalised consensus value
    if (Maxsco>0.0) Consval/=(Maxsco*Pno); else Consval=0.0;

    // get most abundant character
    int Maxno=-1;
    for (j=0; j<Aacodelen; j++)
    {
	if (Anum[j]>Maxno)
	{
	    Maxno=Anum[j]; Aj=Aacodes[j];
	}
	if (Anum[j]==Maxno)	// tie with previous best
	{
	    if (simil(Aj, Aj)<simil(Aacodes[j], Aacodes[j])) Aj=Aacodes[j];    // new is better
	}
    }
    
    delete [] Anum;
    return(Aj);	// the consensus char
}
// END of cons()

/* reset(): brings the object back to the same state it was in when
 * the default ctor created it (26 uppercase letters as code string, 
 * unit matrix as similarity matrix).
 */
void Simil_::reset()
{
    delete [] Aacodes;
    Aacodes=new char [27];
    strcpy(Aacodes, "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    Sim.set_size(26);
    Sim.set_values(); Sim.diag_matrix();    // unit matrix
}
// END of reset()

/* pos(): returns the position of the character C in the Aacodes string.
 * Converts C silently to uppercase. If C is not found then it is changed
 * into 'X' (for unknown) and 'X' is attempted to be found. If unsuccessful, 
 * then a negative number is returned. On success, the index is returned.
 */
int Simil_::pos(char& C) const
{
    C=toupper(C);
    const char *Pos=strchr(Aacodes, C);
    if (Pos!=NULL) return(Pos-Aacodes);	    // success
    
    cerr<<"\n? Simil_::pos('"<<C<<"'): Unknown code, replaced by 'X'\n";
    Pos=strchr(Aacodes, 'X');
    if (Pos!=NULL) return(Pos-Aacodes);	    // success with 'X'
    
    cerr<<"\n? Simil_::pos('X'): Not found\n";
    return(-1);	// error
}
// END of pos()

// ---- Input ----

/* read_file(): reads a similarity matrix and the corresponding amino acid
 * code string from the file Fname. The data members will be changed only
 * upon successful completion. The format is as follows:-
 * Lines beginning with '#' are comments and will be ignored.
 * The first non-comment line will be the amino acid code string. It must
 * contain characters from the set [A-Z\-] only: any other char will trigger
 * an error. The following lines will be interpreted as the rows of the
 * matrix. The dimension is the length of the code string: dim mismatches
 * will generate errors. Note that the format has changed from DRAGON 3.x
 * and the "minimal score item" is no longer supported.
 * Return value: 0 on error,  >0 if OK.
 */
int Simil_::read_file(const char *Fname)
{
    ifstream Infile(Fname);
    if (!Infile)
    {
	cerr<<"\n? Simil_::read_file(\""<<Fname<<"\"): Cannot open\n";
	return(0);
    }
    
    Infile>>(*this);
    return(Infile.good() || Infile.eof()); // done
}
// END of read_file()

/* >>: reads a similarity matrix from Inf. See comments to read_file() */
istream& operator>>(istream& Inf, Simil_& Simil)
{
    if (!Inf)
    {
	cerr<<"\n? >>Simil_: Cannot read from stream\n";
	return(Inf);
    }
    
    char Line[LINELEN];
    memset(Line, 0, LINELEN);
    istrstream Inline(Line, LINELEN-1);
    char *Codes=NULL, *Nlpos;
    Sqmat_ Tempsim;
    int i, ix, jx, Aano, Lineno;
    double Val, Minval=HUGE_VAL;
    
    for (Lineno=0; Inf.good(); Lineno++)
    {
	Line[0]='\0';	// make buffer look empty
	Inf.getline(Line, LINELEN, '\n');
	if (!Inf)
	{
	    cerr<<"\n? >>Simil_: Cannot read line "<<(Lineno+1)<<endl;
	    Inf.clear(Inf.rdstate()|ios::badbit);
	    return(Inf);
	}
	if (Line[0]=='#') continue; // skip comments
	
	if (Codes==NULL)    // haven't seen the code string yet
	{
	    Nlpos=strrchr(Line, '\n');
	    if (Nlpos!=NULL) *Nlpos='\0';   // terminate at \n
	    Aano=strlen(Line);
	    for (i=0; i<Aano; i++)
		if (!isupper(Line[i]))
		{
		    cerr<<"\n? >>Simil_: Non-AA char \'"<<Line[i]
			<<"\' in code string \""<<Line<<"\"\n";
		    Inf.clear(Inf.rdstate()|ios::badbit);
		    return(Inf);
		}
	    
	    // OK, store in Codes
	    Codes=new char [Aano+1];
	    strcpy(Codes, Line);
	    Codes[Aano]='\0';    // paranoia
	    Tempsim.set_size(Aano);	// set temporary similarity matrix size
	    ix=0;   // first row should come
	    continue;	// get new line
	}
	
	// read data lines
	Inline.seekg(ios::beg);
	for (jx=0; jx<Aano; jx++)
	{
	    Inline>>Val;
	    if (!Inline)
	    {
		cerr<<"\n? >>Simil_: Cannot parse line "<<(Lineno+1)<<endl;
		delete [] Codes;
		Inf.clear(Inf.rdstate()|ios::badbit);
		return(Inf);
	    }
	    if (Val<Minval) Minval=Val;	// get minimal value
	    Tempsim[ix][jx]=Val;    // store value temporarily
	}
	
	if (jx<Aano)
	{
	    cerr<<"\n? >>Simil_: Line "<<(Lineno+1)<<" is too short\n";
	    delete [] Codes;
	    Inf.clear(Inf.rdstate()|ios::badbit);
	    return(Inf);
	}
	
	ix++;	// prepare for next line
	if (ix>=Aano) break;	// finish up?
    }
    
    // must be OK if this stage is reached
    delete [] Simil.Aacodes;		// update code string
    Simil.Aacodes=new char [Aano+1];
    strcpy(Simil.Aacodes, Codes);
    Simil.Aacodes[Aano]='\0';
    delete [] Codes;
    
    // shift the whole matrix so that the minimal value is 0.0
    for (ix=0; ix<Aano; ix++)
	for (jx=0; jx<Aano; jx++)
	    Tempsim[ix][jx]-=Minval;
    Simil.Sim=Tempsim;    // copy matrix
    return(Inf);
}
// END of operator>>

ostream& operator<<(ostream& Out, const Simil_& S)
{
    Out<<"Amino acid similarity matrix\n";
    Out<<"Amino acids:"<<S.Aacodes<<endl;
    Out<<"The matrix:\n"<<S.Sim;
    return(Out);
}
// END of operator<<

#undef LINELEN

// ==== END OF METHODS Simil.c++ ====
