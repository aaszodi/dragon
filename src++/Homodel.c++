// ==== PROJECT DRAGON: METHODS Homodel.c++ ====

/* Class for distance-based homology modelling. */

// SGI C++, IRIX 6.2, 30. Apr. 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <fstream.h>
#include <errno.h>
#include <unistd.h> /* for unlink() */
#include <stdio.h>  /* for sprintf() */
#include <float.h>  /* for FLT_MAX */

// ---- MODULE HEADERS ----

#include "Homodel.h"
#ifdef USE_PVM
#include "Pvmtask.h"
#endif

// ---- DEFINITIONS ----

#ifndef FLT_MAX
#define FLT_MAX (1e+20)
#endif

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this: Define NO_MATHFLOATFUNC on the command line
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#endif

// ==== Homodel_ METHODS ====

// ---- Modification ----

/* read_knownstr(): reads the known structure(s) from a PDB file Pdbf.
 * More than one structure may be specified in the file as separate
 * chains. The sequences of them are extracted and compared to the
 * sequences already in the alignment: chains which were not found in
 * the alignment are skipped. The C-alpha coordinates are then stored.
 * Nothing is done and a negative integer is returned on error
 * (inaccessible or unreadable PDB file), otherwise the number of
 * structures successfully identified is returned. 0 is returned
 * if NULL or "" was specified as an input PDB file.
 */
int Homodel_::read_knownstr(const char *Pdbf)
{
    // don't do anything if no file was specified
    if (Pdbf==NULL || !strlen(Pdbf))
    {
	Knownno=0; return(0);
    }
    
    // get the PDB information: all atoms although only the C-alphas are used
    Pdbentry_ *Pdb=get_pdb(Pdbf, ALLATOMS, STRICT);
    if (Pdb==NULL)
    {
	cerr<<"\n? Homodel_::read_knownstr(\""<<Pdbf<<"\"): Cannot read PDB file\n";
	return(-1);
    }
    
    // this PDB file looks OK, make room for the structures
    delete [] Knownstructs; 
    Knownstructs=new Known_ [Pdb->Chainno];
    Bestknown=NULL; Knownno=0;
    
    // process all chains in turn
    unsigned int Chno;
    Chain_ *Chain;
    float Sim=0.0;
    const Align_& Aln=Pol.align();  // get ref to the alignment
    
    for (Chain=Pdb->Chains, Chno=0; Chno<Pdb->Chainno; Chain++, Chno++)
    {
	if (Chain->Type=='X') continue;	// not protein, skip
	
	// find out if the sequence is in the alignment
	char *Alnseq=NULL;	// this is just a buffer
	unsigned int k;
	cout<<"# Sequence of known structure "<<(Chno+1)<<endl;
	for (unsigned int q=0; q<strlen(Chain->Seq); q++)
	{
	    cout<<Chain->Seq[q];
	    if ((q+1) % 60 ==0) cout<<endl; // 60 chars per line
	    else if ((q+1) % 10 ==0) cout<<' ';	// in groups of 10
	}
	cout<<endl;
	for (k=0; k<Aln.seq_no(); k++)
	{
	    Aln.seq(k, Alnseq);
	    if (!strcmp(Chain->Seq, Alnseq)) break;	// OK, matches
	}
	if (k>=Aln.seq_no())
	{
	    delete [] Alnseq;
	    continue;	// struct not found in alignment
	}

	// store this structure
	Knownstructs[Knownno].Seq=Alnseq;   // store sequence
	Knownstructs[Knownno].Seqidx=k;	// store sequence index in alignment
	get_ca(Chain, Knownstructs[Knownno].Cas);   // extract CA coordinates
	Knownstructs[Knownno].Sim=Sim=
	    Pol.seq_simil((int)Pol.master()-1, k);    // similarity to target
	
	// save struct most similar to target so far
	if (Bestknown==NULL || Bestknown->Sim<Sim)
	    Bestknown=Knownstructs+Knownno;
	cout<<"# Chain "<<(Chno+1)<<" of \""<<Pdbf
	    <<"\" is the "<<(k+1)<<". sequence in the alignment\n";
	
	Knownno++;
	
    }	    // for Chain
    
    if (Bestknown!=NULL)
    {
	cout<<"# The "<<(Bestknown->Seqidx+1)
	    <<". sequence in the alignment is the most similar to the target\n";
    }
    else
	cerr<<"\n? Homodel_::read_knownstr(): Cannot find most similar chain\n";
    
    // cleanup
    free_pdb(Pdb);
    return(Knownno);
}
// END of read_knownstr()

/* NOTE: the following method is used for PVM runs only */
#ifdef USE_PVM

/* str_readknown(): reads the PDB structure from the string Pdbstr.
 * This is a horrible kludge: the string is written to a temporary
 * file which is then processed by read_knownstr(). The reason is
 * that the PDB reader is in C, and it wouldn't be easy to convert
 * all those scanf()-s to sscanf()-s. This method is provided
 * for the PVM inter-process message passing.
 * Return values are exactly the same as in read_knownstr().
 */
int Homodel_::str_readknown(const char *Pdbstr)
{
    // don't do anything if Pdbstr is empty
    if (Pdbstr==NULL || !strlen(Pdbstr))
    {
	cerr<<"\n? Homodel_::str_readknown(): PDB string is empty\n";
	Knownno=0; return(0);
    }
    
    extern Pvmtask_ Pvmtask;
    if (!Pvmtask.is_slave())
    {
	cerr<<"\n? Homodel_::str_readknown(): Not running as a PVM slave\n";
	return(-4);
    }
    
    /* Set up a name for the tempfile: the PVM task ID string is used
     * for this purpose (which must be unique). 
     */
    char Pdbtemp[256];
    ofstream Tmpf;
    
    sprintf(Pdbtemp, "D4pdbtemp_%s", Pvmtask.id_str());
    cout<<"# Homodel tempfile name:"<<Pdbtemp<<endl;
    Tmpf.open(Pdbtemp);
    if (!Tmpf)
    {
	cerr<<"\n? Homodel_::str_readknown(): Tempfile \""
	    <<Pdbtemp<<"\" cannot be opened\n";
	return(-5);
    }
    
    // Now comes the horror: write Pdbstr to tempfile and read back
    Tmpf.write(Pdbstr, strlen(Pdbstr));	// unformatted write
    Tmpf.close();
    int Retval=read_knownstr(Pdbtemp);
    
    if (unlink(Pdbtemp)<0)  // remove tempfile
    {
	cerr<<"\n? Homodel_::str_readknown(): Removal of tempfile \""
	    <<Pdbtemp<<"\" failed. Never mind\n";
    }
    return(Retval); // read_knownstr()'s result
}
// END of str_readknown()
#endif	/* USE_PVM */

/* get_ca(): updates the point set Cas
 * which holds the coordinates of the CA (C-alpha) atoms. Private static
 */
void Homodel_::get_ca(const Chain_ *Chain, Points_& Cas)
{
    // adjust size and set full activation
    Cas.len_dim(Chain->Aano, 3);
    
    Atom_ *Cur;
    int i, j;
    
    /* trundle along the chain: the C-alpha always comes
     * before the corresponding side chain (we hope)
     */
    for (Cur=Chain->Atoms, i=j=0; i<Chain->Atomno; i++, Cur++)
    {
	// skip all atoms except C-alphas
	if (strcmp(Cur->Id, "CA") || !(Cur->Alt==' ' || Cur->Alt=='A'))
	    continue;
	
	// store C-alpha coords
	Cas[j][0]=Cur->X; Cas[j][1]=Cur->Y; Cas[j][2]=Cur->Z; 
	j++; 
    }
}
// END of get_ca()
    
// ---- Restraint generation ----

/* make_restrs(): builds a list of CA distance restraints
 * for all residue pairs in the known structure(s) which participate
 * in the alignment and are closer than the Maxdist (not squared).
 * The lower bound is the shortest distance found
 * in the known structures, the upper bound is the largest.
 * The residue pairs must be separated by Minsepar residues in
 * the sequence (Minsepar==2 is the minimum, allowing for
 * (i, i+2) pairs). 2 is also the default.
 * From V4.11.2 on, also prepares for RMS checks between model
 * structure and the best known structure (to get the right enantiomer).
 */
List1_<Restr_> Homodel_::make_restrs(float Maxdist, int Minsepar)
{
    List1_<Restr_> Rlist;
    if (!Knownno)
    {
	cerr<<"\n? Homodel_::make_restrs(): No known structure\n";
	return(Rlist);	// returns empty list
    }
    
    /* Prepare for the RMS comparison between the best known structure
     * and the model. While looking for the close pairs, 
     * equivalence masks are also constructed. Two residues in
     * the best scaffold and the model are considered equivalent if
     * they are in the same alignment position. The comparison
     * weights come from the conservation values at these positions.
     */
    
    /* Set weight array and mask lengths. Note that 
     * the model has 2 more "pseudo-C-alphas" at both
     * ends of the chain (NH3+, COO-)
     */
    unsigned int Len=Pol.len();	// master seq length
    Weight.dim(Len);
    Knownmask.len(Bestknown->Cas.len());    // use most homologous
    Knownmask.set_values(false);
    Modelmask.len(Len+2);   // 2 added for N/C termini
    Modelmask.set_values(false);
    
    register unsigned int k, w;
    register int ai, aj, mi, mj; // alignment and master idx
    int Bestidx=Bestknown-Knownstructs;	// index of best known structure
    register float Ci, Cj, D2, D2low, D2hi;
    Restr_ R;
    const Align_& Aln=Pol.align();
    int *Si=new int [Knownno];	// seq position array for known structs
    int *Sj=new int [Knownno];
    
    // use the squared maximal distance
    cout<<"# Maximal distance for homology-derived restraints: "<<Maxdist<<endl;
    Maxdist*=Maxdist;
    
    cout<<"# List of equivalent residues\n\n# TARG";
    for (k=0; k<Knownno; k++)
    {
	cout<<'\t'<<(Knownstructs[k].Seqidx+1);
	if (k==Bestidx) cout<<'*';
    }
    cout<<endl;
    for (k=0; k<=Knownno; k++) cout<<"--------";
    cout<<endl;
    
    // set minimal sequential separation
    if (Minsepar<2)
    {
	cerr<<"\n? Homodel_::make_restrs(): Minsepar="<<Minsepar<<"<2, set to 2\n";
	Minsepar=2;
    }
    
    // count restraints here
    const unsigned int SEPAR_CATEGORIES=5;
    enum Separ_ { VERY_CLOSE=4, CLOSE=10, MEDIUM=20, DISTANT=50};
    unsigned int *Separs=new unsigned int [5];
    int Sep;
    for (k=0; k<SEPAR_CATEGORIES; k++) Separs[k]=0;
    
    // scan all residue pairs
    for (mi=w=0; mi<Len; mi++)
    {
	/* The mi:th position of the master sequence is in
	 * the ai:th position of the alignment. Si[k] will hold
	 * the squence position of the k:th known structure
	 * corresponding to the ai:th alignment column or -1
	 * if the k:th struct is gapped there. Easy :-)
	 */
	// ai==mi if the master is the consensus
	ai=(!Pol.master())? mi: Aln.align_pos(Pol.master()-1, mi);
	if (ai<0) continue; // error
	Ci=Pol.cons(mi);    // conservation value
	
	// get the residue type and index for each known struct in ai:th align pos: <0 if gap
	cout<<Pol.aa(mi)<<"["<<setw(3)<<(mi+1)<<"]";
	for (k=0; k<Knownno; k++)
	{
	    Si[k]=Aln.seq_pos(Knownstructs[k].Seqidx, ai);
	    if (Si[k]>=0) cout<<'\t'<<Knownstructs[k].Seq[Si[k]]<<"["<<setw(3)<<(Si[k]+1)<<"]";
	    else cout<<"\t------";
	}
	cout<<endl;
	
	if (Si[Bestidx]>=0)
	{
	    // found equiv. res. in best known for mi:th in master
	    Knownmask.set_bit(Si[Bestidx]);
	    Modelmask.set_bit(mi+1);    // shift: [0]th pos is NH3+
	    Weight[w++]=Ci;	// comparison weight for Si[..]:[mi] pair
	}
	
	// same for the other half of the restraint pairs
	for (mj=mi+Minsepar; mj<Len; mj++)
	{
	    aj=(!Pol.master())? mj: Aln.align_pos(Pol.master()-1, mj);
	    if (aj<0) continue; // error
	    Cj=Pol.cons(mj);    // conservation value
	    
	    // get the residue index for each known struct in aj:th align pos: <0 if gap
	    for (k=0; k<Knownno; k++)
		Sj[k]=Aln.seq_pos(Knownstructs[k].Seqidx, aj);
	    
	    // construct the restraint between the mi:th and mj:th residues
	    R.pos(1, mi+1); R.pos(2, mj+1);	// 1..Rno indexing
	    R.atom(1, "CA"); R.atom(2, "CA");	// CA:CA only
	    R.strict(sqrtf(Ci*Cj)); // strictness is geom.avg. of conservation
	    
	    // find shortest and longest observed distances in the known structs
	    D2low=FLT_MAX; D2hi=-D2low;
	    for (k=0; k<Knownno; k++)
	    {
		if (Si[k]<0 || Sj[k]<0)
		    continue;	// k:th had a gap either w/ mi or mj, no good
		
		// get actual squared CA:CA distance from k-th known 
		D2=diff_len2(Knownstructs[k].Cas[Si[k]], Knownstructs[k].Cas[Sj[k]]);
		if (D2>Maxdist)
		    continue;	// not close enough
		
		if (D2<D2low) D2low=D2;
		if (D2>D2hi) D2hi=D2;
	    }
	    
	    if (D2low==FLT_MAX || D2hi==-FLT_MAX)
		continue;   // all known structs had gaps here or not close enough
	    
	    // widen the range by +/- 5 %
	    D2low*=0.9025; D2hi*=1.1025;
	    
	    R.low2(D2low); R.up2(D2hi);	// store restraints at last
	    Rlist+=R;
	    
	    // update separation categories (not terribly elegant)
	    Sep=abs(int(mj)-int(mi));
	    if (Sep<=VERY_CLOSE) Separs[0]++;
	    else if (Sep<=CLOSE) Separs[1]++;
	    else if (Sep<=MEDIUM) Separs[2]++;
	    else if (Sep<=DISTANT) Separs[3]++;
	    else Separs[4]++;
	    
	}	// for mj
    }	    // for mi
    
    // list separation "statistics"
    cout<<"# Restraint distribution by sequential separation\n";
    cout<<"2.."<<VERY_CLOSE<<":\t"<<Separs[0]<<endl;
    cout<<(VERY_CLOSE+1)<<".."<<CLOSE<<":\t"<<Separs[1]<<endl;
    cout<<(CLOSE+1)<<".."<<MEDIUM<<":\t"<<Separs[2]<<endl;
    cout<<(MEDIUM+1)<<".."<<DISTANT<<":\t"<<Separs[3]<<endl;
    cout<<DISTANT<<"+:\t"<<Separs[4]<<endl;
    
    Weight.dim(w);  // readjust weight vector's length to no. of equivalent pos.
    delete [] Separs; delete [] Si; delete [] Sj;
    return(Rlist);
}
// END of make_restrs()

/* hand_check(): compares the model C-alpha coordinates in Model to
 * the C-alpha coordinates of the scaffold structure most homologous to the model.
 * Returns 1 if Model is more similar to the scaffold than its mirror image, 
 * -1 if a flip was needed (which is done inside) and 0 on error
 * or if Model was not 3D. Cf. the hand_check() fn in "Sterchem".
 */
int Homodel_::hand_check(Points_& Model)
{
    if (Model.dim()!=3 || Bestknown==NULL) 
	return(0);  // not in 3D or not in homology mode, do nothing

    Points_ &Cas=Bestknown->Cas;    // coords of the most homologous structure
    
    // set the equivalence masks
    Bits_ Oldmodelmask=Model.mask(Modelmask);
    Cas.mask(Knownmask); 
    
    // center the masked structures (inactive points stay the same!!!)
    Vector_ Modctr=Model.centroid(Weight); Model-=Modctr;
    Vector_ Casctr=Cas.centroid(Weight); Cas-=Casctr;
    
    // do the unflipped match
    if (!Hr.best_rot(Cas, Model, Weight))
    {
	cerr<<"\n? Homodel_::hand_check(): Rank deficiency in unflipped rotation\n";
	Model+=Modctr; Model.mask(Oldmodelmask);
	Cas+=Casctr; Cas.mask(true);
	return(0);
    }
    double Rms=Hr.get_rms(Cas, Model, Weight);
    if (Rms<0.0)
    {
	cerr<<"\n? Homodel_::hand_check(): Cannot get RMS in unflipped rotation\n";
	Model+=Modctr; Model.mask(Oldmodelmask);
	Cas+=Casctr; Cas.mask(true);
	return(0);
    }
    
    // now do the flipped match
    Model+=Modctr; Model.mask(true);	// switch all points on: Model moved back
    Points_ Flipmodel(Model);
    Flipmodel*=-1.0;    // invert the model (all coords)
    Flipmodel.mask(Modelmask);	// apply equivalence mask
    Modctr=Flipmodel.centroid(Weight); Flipmodel-=Modctr; // centre active flipped
    
    Model.mask(Oldmodelmask);	// on error, the original Model is retained
    
    if (!Hr.best_rot(Cas, Flipmodel, Weight))
    {
	cerr<<"\n? Homodel_::hand_check(): Rank deficiency in flipped rotation\n";
	Cas+=Casctr; Cas.mask(true); return(0);
    }
    double Rmsflip=Hr.get_rms(Cas, Flipmodel, Weight);
    if (Rmsflip<0.0)
    {
	cerr<<"\n? Homodel_::hand_check(): Cannot get RMS in flipped rotation\n";
	Cas+=Casctr; Cas.mask(true); return(0);
    }
    Cas+=Casctr; Cas.mask(true);    // reset Calphas of known
    
    // choose the enantiomer with the lower RMS value
    cout<<"HAND: (homol) RMS="<<Rms<<", FLIP="<<Rmsflip<<endl;
    if (Rms<=Rmsflip)
    {
	// keep the original
	return(1);
    }
    else
    {
	// use the flipped enantiomer
	Flipmodel+=Modctr; Flipmodel.mask(true);
	Model=Flipmodel; Model.mask(Oldmodelmask);
	return(-1);
    }
}
// END of hand_check()

// ==== END OF METHODS Homodel.c++ ====
