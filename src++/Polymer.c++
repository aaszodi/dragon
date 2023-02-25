// ==== PROJECT DRAGON: METHODS Polymer.c++ ====

/* The Polymer_ class holds all information about the model
 * chain,  viz. sequence, hydrophobicity, conservation, 
 * accessibility etc.*/

// SGI C++, IRIX 6.2, 2. Sept. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <iomanip.h>
#include <strstream.h>

// ---- CLASS HEADER ----

#include "Polymer.h"

// ---- UTILITY HEADER ----

#include "Stat2.h"

// ---- DEFINITIONS ----

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#ifndef CA_BUMP
#define CA_BUMP (2.0)	/* C-alpha bump radius (not squared) */
#endif

// ==== Polymer_ METHODS ====

// ---- Static initialisation ----

Property_ Polymer_::Hyphob(Property_::Hyphobdef);
Property_ Polymer_::Volume(Property_::Volumedef);
Simil_ Polymer_::Simil;
Distpred_ Polymer_::Dp;
Acdist_ Polymer_::Acdist;

// ---- Access ----

/* estim_dist(): returns the NON-SQUARED estimated distance between
 * residues R1, R2 estimated from their conserved hydrophobicity
 * scores. Prints a warning and returns 0.0 if R1, R2 are out of
 * range [0..Rno-1].
 */
double Polymer_::estim_dist(unsigned int R1, unsigned int R2)
{
    // shall we do the estimation?
    if (Changed)
    {
	for (unsigned int i=0; i<len(); i++)
	    Consphob[i]=Monomers[i].Cons*Monomers[i].Phob;
	Dp.estim_params(Consphob);
	Changed=0;
    }
    if (!len() || R1>len() || R2>len())
    {
	cerr<<"\n? Polymer_::estim_dist("<<R1<<", "<<R2<<"): Invalid index, 0.0 returned\n";
	return(0.0);
    }
    return(Dp.dist_phob(Consphob[R1]+Consphob[R2]));
}
// END of estim_dist()

/* master(Mseq): changes the master sequence within the alignment. An
 * argument of 0 means the consensus, otherwise the Mseq-1:th sequence
 * will be the consensus. No action is taken if Mseq is invalid.
 * Returns old master sequence number.
 */
unsigned int Polymer_::master(unsigned int Mseq)
{
    if (!Align.seq_no())
    {
	cerr<<"\n? Polymer_::master("<<Mseq<<"): No sequences, no action\n";
	return(0);
    }
    if (Mseq>=Align.seq_no())
    {
	cerr<<"\n? Polymer_::master("<<Mseq
	    <<"): Invalid master sequence number, consensus used\n";
	Mseq=0;
    }
    if (Master==Mseq) return(Master);	// no change, do nothing
    
    unsigned int Oldm=Master; Master=Mseq;
    update_members(ALIGN);  // equivalent to alignment change
    return(Oldm);
}
// END of master()

/* seq_simil(): calculates a rough sequence similarity between
 * the S1:th and S2:th sequences in the current alignment using
 * the current similarity matrix. S1, S2 in [0..Seqno]. If S1 or S2
 * are -1, then the consensus sequence of the alignment is used.
 * The score returned is an average of pairwise scores between all non-gapped
 * positions. Returns -1.0 on error.
 */
float Polymer_::seq_simil(int S1, int S2) const
{
    if (S1<-1 || S1>=Align.seq_no() || S2<-1 || S2>=Align.seq_no())
    {
	cerr<<"\n? Polymer_::seq_simil("<<S1<<","<<S2<<"): Invalid sequence(s)\n";
	return(-1.0);
    }
    
    const char *Posstr;
    unsigned int k, n=0;
    float Sim=0.0;
    char A, B;
    double Consval; // dummy consensus value required by Simil_::cons()
    
    // walk along the alignment
    for (k=0; k<Align.len(); k++)
    {
	Posstr=Align.pos(k);	// get the k-th alignment position
	
	// obtain the amino acid codes (S1,S2==-1: then consensus code)
	A=(S1==-1)? Simil.cons(Posstr, Consval): Posstr[S1];
	if (A=='-') continue;	// skip gaps
	B=(S2==-1)? Simil.cons(Posstr, Consval): Posstr[S2];
	if (B=='-') continue;

	Sim+=Simil.simil(A, B);	// sum similarity values
	n++;	// count matches
    }
    
    if (!n)
    {
	cerr<<"\n? Polymer_::seq_simil("<<S1<<","<<S2<<"): No matches\n";
	return(0.0);
    }
    else return(Sim/n);
}
// END of seq_simil()

// ---- Input ----

/* read_aln(): attempts to read an alignment file from Fname. If Mseq==0,
 * then the "master sequence" will be the consensus of the alignment (default);
 * otherwise, the (Mseq-1)-th sequence from the alignment will be the
 * master sequence (i.e. the sequence of the model chain). If the input
 * is not successful or there are not enough sequences then nothing will
 * be changed; however, a successful alignment input will cause all the
 * other fields be updated.
 * Return value: 0 on failure, the new sequence length on success.
 */
int Polymer_::read_aln(const char *Fname, unsigned int Mseq)
{
    if (!Align.read_file(Fname)) return(0);
    
    // check choice of master sequence
    if (!Mseq && Mseq>=Align.seq_no())
    {
	cerr<<"\n? Polymer_::read_aln("<<Fname<<","<<Mseq
	    <<"): Invalid master sequence number, consensus used\n";
	Mseq=0;
    }
    Master=Mseq;
    
    update_members(ALIGN);
    return(Monomers.len());
}
// END of read_aln()

/* str_aln(): does the same thing as read_aln() but reads from 
 * a string Str.
 */
int Polymer_::str_aln(const char *Str, unsigned int Mseq)
{
    istrstream Ifs(Str); Ifs>>Align;
    if (!Align.len()) return(0);
    
    // check choice of master sequence
    if (!Mseq && Mseq>=Align.seq_no())
    {
	cerr<<"\n? Polymer_::str_aln(...,"<<Mseq
	    <<"): Invalid master sequence number, consensus used\n";
	Mseq=0;
    }
    Master=Mseq;
    
    update_members(ALIGN);
    return(Monomers.len());
}
// END of read_aln()

/* read_phob(), read_vol(), read_acdist, read_simil(): 
 * read new data from a file Fname. Update the
 * sequence records accordingly. Return 0 on failure, 1 on success.
 * The corresponding str_* methods do exactly the same but
 * read the input from a string Str rather than from a file.
 */
int Polymer_::read_phob(const char *Fname)
{
    if (!Hyphob.read_file(Fname)) return(0);	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(HYPHOB);
    return(Monomers.len());
}
// END of read_phob()

int Polymer_::str_phob(const char *Str)
{
    istrstream Ifs(Str); Ifs>>Hyphob;
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(HYPHOB);
    return(Monomers.len());
}
// END of str_phob()

int Polymer_::read_vol(const char *Fname)
{
    if (!Volume.read_file(Fname)) return(0);	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(VOLUME);
    return(Monomers.len());
}
// END of read_vol()

int Polymer_::str_vol(const char *Str)
{
    istrstream Ifs(Str); Ifs>>Volume;	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(VOLUME);
    return(Monomers.len());
}
// END of str_vol()

int Polymer_::read_acdist(const char *Fname)
{
    if (!Acdist.read_file(Fname)) return(0);	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(ACDIST);
    return(Monomers.len());
}
// END of read_acdist()

int Polymer_::str_acdist(const char *Str)
{
    istrstream Ifs(Str); Ifs>>Acdist;	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(ACDIST);
    return(Monomers.len());
}
// END of str_acdist()

int Polymer_::read_simil(const char *Fname)
{
    if (!Simil.read_file(Fname)) return(0);	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(SIMIL);
    return(Monomers.len());
}
// END of read_simil()

int Polymer_::str_simil(const char *Str)
{
    istrstream Ifs(Str); Ifs>>Simil;	// input
    if (!Align.seq_no()) return(1);  // no sequences to be modified

    update_members(SIMIL);
    return(Monomers.len());
}
// END of str_simil()

/* update_members(): called after a successful disk file read which
 * modified one of the members. The modified member is indicated by
 * the Mod parameter. Depending on the nature of this, some other members
 * must be updated (such as the corresponding fields in the elements of
 * the Monomers array). Private
 */
void Polymer_::update_members(int Mod)
{
    if (!Mod) return;   // no action
    
    const char *Posstr;
    register int i, k, Rno;
    double Ftemp, Cons;
    Stat_ Cstat;
    
    if (Master)	    // use one of the sequences as polymer sequence
    {
	Rno=Align.seq_len(Master-1);
	char Aa;
	double Ftemp;
	
	if (Mod & ALIGN) Monomers.len(Rno);  // new sequence, realloc
	
	for (i=k=0; i<Align.len(); i++, k++)
	{
	    Posstr=Align.pos(i);
	    if ((Aa=Posstr[Master-1])=='-') { k--; continue; }	// skip gaps
	    
	    if (Mod & ALIGN)
		Monomers[k].Aa=Aa;	// amino acid type
	    if (Mod & (ALIGN|SIMIL))
	    {
		Simil.cons(Posstr, Cons);   // consensus
		Monomers[k].Cons=Cons;
		Cstat+=Cons;
	    }
	    if (Mod & (ALIGN|HYPHOB))
		Monomers[k].Phob=Hyphob[Aa];	// individual phobicity
	    if (Mod & (ALIGN|VOLUME))
	    {
		Monomers[k].Bumpb=Ftemp=
		    pow(3.0*Volume[Aa]/(4.0*M_PI), 1.0/3.0); // C-beta bump radius
		Ftemp+=CA_BUMP;	    // add C-alpha bump radius (bump dist)
		Monomers[k].Bumpab=Ftemp*Ftemp;	// store squared AB bump dist
	    }
	    if (Mod & (ALIGN|ACDIST))
	    {
		Ftemp=Acdist.scc_dist(Aa, "CA");
		Monomers[k].Abdist=Ftemp*Ftemp;	// CA:SCC dist squared
	    }
	}
    }
    else    // use consensus
    {
	Rno=Align.len();
	if (Mod & ALIGN) Monomers.len(Rno);  // new sequence, realloc
	
	for (i=0; i<Rno; i++)
	{
	    Posstr=Align.pos(i);
	    if (Mod & (ALIGN|SIMIL))
	    {
		Monomers[i].Aa=Simil.cons(Posstr, Cons);
		Monomers[i].Cons=Cons;
		Cstat+=Cons;
	    }
	    if (Mod & (ALIGN|SIMIL|HYPHOB))
		Monomers[i].Phob=Hyphob.avg_val(Posstr);	// average phobicity
	    if (Mod & (ALIGN|SIMIL|VOLUME))
	    {
		Monomers[i].Bumpb=Ftemp=
		    pow(3.0*Volume.avg_val(Posstr)/(4.0*M_PI), 1.0/3.0);
		Ftemp+=CA_BUMP;
		Monomers[i].Bumpab=Ftemp*Ftemp;
	    }
	    if (Mod & (ALIGN|SIMIL|ACDIST))
	    {
		Ftemp=Acdist.scc_dist(Monomers[i].Aa, "CA");
		Monomers[i].Abdist=Ftemp*Ftemp;	// CA:SCC average dist squared
	    }
	}
    }
    
    // get the average and SD of the conservation
    if (Mod & (ALIGN|SIMIL))
    {
	Cavg=Cstat.avg(); Csd=Cstat.sd();
    }
    
    if (Mod & ALIGN) Consphob.len(Rno);	// realloc if necessary
    
    /* NOTE: the only update that is deferred is the nonlinear
     * parameter estimation for the phobicity->distance prediction.
     * A flag is set here and will be processed by estim_dist().
     */
    if (Mod & (ALIGN|SIMIL|HYPHOB))
	Changed=1;
}
// END of update_members()

// ---- Output ----

ostream& operator<<(ostream& Out, const Polymer_& P)
{
    Out<<"# No. of sequences = "<<P.Align.seq_no()<<", model = ";
    if(P.Master) Out<<"Seq. #"<<P.Master;
    else Out<<"consensus";
    Out<<", no. of residues = "<<P.len()<<endl;
    if (!P.len()) return(Out);
    Out<<P.Align<<endl;
    
    Out<<"#\tTarget\tCons\tPhob\tBrad\tAcdist\tAlignment to target\n";
    Out.precision(3);
    int k;
    for (unsigned int i=0; i<P.len(); i++)
    {
        if (P.Master) k=P.align().align_pos(P.Master-1,i);
	Out<<(i+1)<<'\t'<<P.aa(i)<<'\t'<<P.cons(i)<<'\t'
	    <<P.phob(i)<<'\t'<<P.bumpb(i)<<'\t'<<sqrt(P.abdist(i))<<'\t'
	    <<P.align().pos(P.Master? k: i)<<endl;
    }
    Out<<"# Average conservation="<<P.cons_avg()<<", SD="<<P.cons_sd()<<endl;
    return(Out);
}

// ==== END OF METHODS Polymer.c++ ====
