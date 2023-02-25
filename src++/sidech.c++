// ==== PROJECT DRAGON: PROGRAM sidech.c++ ====

/* Decorates a protein main chain with side chains,
 * given a MULTAL alignment, a main chain (in PDB format) and
 * another PDB file containing homologous structures.
 * The main chain comes from the "catomain" program
 * that adds peptide bonds to C-alpha chains.
 */

// SGI C++ 4.0, IRIX 5.3, 30. Apr. 1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

// ---- MODULE HEADERS ----

#include "Aacid.h"
#include "Array.h"
#include "Align.h"
#include "Hirot.h"
#include "pdbprot.h"

// ---- CLASSES ----

/* struct Hom_: stores the AA description of a chain homologous
 * to the target, together with the sequence index in the alignment.
 */
struct Hom_
{
    int Seqidx;	// position in alignment, -1 if not there
    Array_<Aacid_> Aas;	// amino acid descriptions (length 0 if not there)
    
    Hom_(): Seqidx(-1), Aas(0) {}
};

// ---- PROTOTYPES ----

static int get_seqpos(const Align_& Align, const Chain_* Chain);
static int get_homs(const Align_& Align, const Pdbentry_* Entry, Array_<Hom_>& Homs);
static int get_aas(const Chain_ *Chain, Array_<Aacid_>& Aas);

static void make_sidechains(const Align_& Align, Hom_& Target, 
	Array_<Hom_>& Homs);
static void equiv_atoms(const Aacid_& Hom, Aacid_& Target);

static void pdb_list(const char *Pdbfname, const char *Alignname, 
	const char *Modelname, const char *Homname, 
	const Array_<Aacid_>& Model);

// ==== MAIN ====

/* The program takes the following parameters: a vertical
 * MULTAL(-like) alignment name, a PDB filename containing
 * a complete N/CA/C/O peptide chain (sidechains are ignored)
 * and another PDB file containing chains homologous to the
 * main chain. Output is to the fourth parameter.
 */
int main(int argc, char *argv[])
{
    if (argc<5)
    {
	cerr<<"\n! Usage: "<<argv[0]<<" alignment mainchain homstruct outfile\n";
	exit(EXIT_FAILURE);
    }
    
    // get the alignment
    Align_ Align;
    if (!Align.read_file(argv[1]))
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read alignment file \""<<argv[1]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    
    // get the main chain, check if it is in the alignment
    Pdbentry_ *Main=NULL;
    if (NULL==(Main=get_pdb(argv[2], ALLATOMS, STRICT)))
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read main chain PDB file \""<<argv[2]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    Hom_ Target;
    Target.Seqidx=get_seqpos(Align, Main->Chains); // first chain only
    if (Target.Seqidx<0)
    {
	cerr<<"\n! "<<argv[0]<<": Main chain PDB file \""<<argv[2]<<"\" is not in alignment\n";
	exit(EXIT_FAILURE);
    }
    get_aas(Main->Chains, Target.Aas);
    free_pdb(Main);
    
    // read the homologous structures
    Pdbentry_ *Homstruct=NULL;
    if (NULL==(Homstruct=get_pdb(argv[3], ALLATOMS, STRICT)))
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read homolog PDB file \""<<argv[3]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    Array_<Hom_> Homs;
    if (!get_homs(Align, Homstruct, Homs))
    {
	cerr<<"\n! "<<argv[0]<<": None of the chains in PDB file \""<<argv[3]<<"\" are in alignment\n";
	exit(EXIT_FAILURE);
    }
    free_pdb(Homstruct);
    
    // perform the side chain fitting
    make_sidechains(Align, Target, Homs);
    
    // output
    pdb_list(argv[4], argv[1], argv[2], argv[3], Target.Aas);
    exit(EXIT_SUCCESS);
}

// ==== FUNCTIONS ====

// ---- Input and setup ----

/* get_seqpos(): given an alignment Align and a PDB chain Chain,
 * -1 is returned if the sequence of Chain is not found in Align, 
 * an index between 0 and Align.seq_no()-1 if found.
 */
static int get_seqpos(const Align_& Align, const Chain_* Chain)
{
    if (Chain->Type=='X') return(-1);	// not protein
    
    // find out if the sequence is in the alignment
    char *Alnseq=NULL;	// this is just a buffer
    int k;
    for (k=0; k<Align.seq_no(); k++)
    {
	Align.seq(k, Alnseq);
	if (!strcmp(Chain->Seq, Alnseq)) break;	// OK, matches
    }
    delete [] Alnseq;
    return((k>=Align.seq_no())? -1: k);
}
// END of get_seqpos()

/* get_homs(): try to find the sequence of the chains in Entry
 * in Align: put the index and AA description for each of
 * those found into Homs.
 * Return value: the number of sequences found (Homs.len()).
 */
static int get_homs(const Align_& Align, const Pdbentry_* Entry, Array_<Hom_>& Homs)
{
    Homs.len(Entry->Chainno);
    
    int k, Idx, Found=0;
    for (k=0; k<Entry->Chainno; k++)
    {
	Idx=get_seqpos(Align, Entry->Chains+k);
	if (Idx>=0)
	{
	    Homs[Found].Seqidx=Idx;
	    get_aas(Entry->Chains+k, Homs[Found].Aas);
	    ++Found;
	}
    }
    Homs.len(Found);
    return(Found);
}
// END of get_homs()

/* get_aas(): converts the PDB chain Chain to the
 * amino acid description array Aas.
 * Return value: the number of residues found, 0 on error.
 */
static int get_aas(const Chain_ *Chain, Array_<Aacid_>& Aas)
{
    if (Chain->Type=='X' || !Chain->Aano) return(0);
    Aas.len(Chain->Aano);
    
    Atom_ *Atoms=Chain->Atoms;
    int Atomidx, Aaidx=-1, Prevresno=-9999;
    Vector_ *Coptr=NULL;
    char Prevrid='~';
    for (Atomidx=0; Atomidx<Chain->Atomno; Atomidx++)
    {
	if (Atoms[Atomidx].Resno!=Prevresno || Atoms[Atomidx].Rid!=Prevrid)
	{
	    // start of new amino acid
	    ++Aaidx;
	    Prevresno=Atoms[Atomidx].Resno;
	    Prevrid=Atoms[Atomidx].Rid;
	    Aas[Aaidx].res_id(Atoms[Atomidx].Aa);   // set amino acid type
	    Aas[Aaidx].mask(false); // does not contain anyone yet
	}
	
	// switch on current atom if possible
	if (!Aas[Aaidx].active(Atoms[Atomidx].Id, true))
	{
	    cerr<<"\n? get_aas(): There's no atom called \""
		<<Atoms[Atomidx].Id<<"\" in amino acid \'"
		<<Aas[Aaidx].res_id()<<"\', skipped\n";
	    continue;
	}
	
	// copy coordinates
	Coptr=Aas[Aaidx].atom(Atoms[Atomidx].Id);
	Coptr->operator[](0)=Atoms[Atomidx].X;
	Coptr->operator[](1)=Atoms[Atomidx].Y;
	Coptr->operator[](2)=Atoms[Atomidx].Z;
    }
    return(Aas.len());
}
// END of get_aas()

// ---- Side chain modelling ----

/* make_sidechains(): establishes equivalences between amino acids
 * in Target and the homologous structures stored in the Homs array, 
 * using the alignment Align. Equivalent atoms for each equivalent
 * residue are identified and the common atom subset in the homologous
 * residues is aligned so that their main-chain portion (N-CA-C-O)
 * overlaps with the target main chain. The averaged sidechain
 * coordinates are then transferred to the target.
 */
static void make_sidechains(const Align_& Align, Hom_& Target, 
	Array_<Hom_>& Homs)
{
    int ti, ai, hi, k, Eqno;
    Bits_ Tmask, Hmask, Ctmask;
    Array_<int> Equiv(Homs.len());
    Vector_ Tctr(3), Hctr(3);	// centroids
    Vector_ W(4);   // main-chain portion rotation weight
    Hirot_ Hr;
    
    /* prepare rotation weight: we happen to know that
     * the order is N, CA, C, O. CA has weight 1, N and C 0.5, 
     * O 0.2 (Willie says these are less reliable)
     */
    W[0]=W[2]=0.5; W[1]=1.0; W[3]=0.2;
    
    // walk along each amino acid in the target
    for (ti=0; ti<Target.Aas.len(); ti++)
    {
	Aacid_& Taa=Target.Aas[ti];

	// if mainchain is not complete, don't bother
	if (!Taa.main_chain())
	{
	    cerr<<"\n? make_sidechains(): No full main chain in target amino acid "
		<<Taa.res_id()<<"-"<<(ti+1)<<endl;
	    continue;
	}
	
	/* The ti:th position of the target sequence is in
	 * the ai:th position of the alignment.
	 */
	ai=Align.align_pos(Target.Seqidx, ti);
	if (ai<0) continue; // error (shouldn't happen)
	Tmask=Taa.mask(true);	// switch all atoms ON
	
	// get equivalent residues if any
	Equiv.set_values(-1);
	for (hi=Eqno=0; hi<Homs.len(); hi++)
	{
	    Equiv[hi]=Align.seq_pos(Homs[hi].Seqidx, ai);
	    if (Equiv[hi]>=0) ++Eqno;	// count non-gapped positions
	}
	
	if (!Eqno)
	{
	    // no equivalent AAs (all gapped)
	    Taa.mask(Tmask);   // reset old activation mask
	    continue;
	}
	
	// generate the overall mask of equivalent atoms in Taa
	for (hi=0; hi<Homs.len(); hi++)
	{
	    if (Equiv[hi]<0) continue;
	    
	    Aacid_& Haa=Homs[hi].Aas[Equiv[hi]];
	    if (!Haa.main_chain())
	    {
		cerr<<"\n? make_sidechains(): No full main chain in homologous structure "
		    <<hi<<", "<<Haa.res_id()<<"-"<<(Equiv[hi]+1)<<endl;
		continue;
	    }
	    equiv_atoms(Haa, Taa);
	}
	
	// rotation
	Ctmask=Taa.mask();	// common mask on target
	Taa.side_chain(false);   // only main-chain active
	Tctr=Taa.centroid(W);	    // center on main-chain centroid
	Taa-=Tctr;
	Taa.mask(Ctmask);	// switch back to common mask
	
	for (hi=0; hi<Homs.len(); hi++)
	{
	    if (Equiv[hi]<0) continue;
	    
	    Aacid_& Haa=Homs[hi].Aas[Equiv[hi]];
	    
	    // mask all equiv residues back with the common mask in target
	    equiv_atoms(Taa, Haa);
	    Hmask=Haa.mask();   // version of Ctmask
if (Ctmask.on_no()!=Hmask.on_no())
{
    cerr<<"Target["<<ti<<"]:\n"<<Taa<<Ctmask;
    cerr<<"Hom["<<hi<<"]:\n"<<Haa<<Hmask;
}

	    // switch off side chains for rotation (main-chain used only)
	    Taa.side_chain(false);
	    Haa.side_chain(false);
	    
	    Hctr=Haa.centroid(W);	// center
	    Haa-=Hctr;
	    
	    // rotate homologous main-chain onto target
	    Hr.best_rot(Haa, Taa, W);
	    Taa.mask(Ctmask);	// switch target to common mask
	    
	    Haa+=Hctr;	// main-chain part back
	    Haa.mask(Hmask);	// main+side ON
	    Haa-=Hctr;	// center main+side to prev. position
	    Haa*=Hr.rot_matrix();
	}
	
	/* Here we have all homologous residues corresponding to the
	 * ti:th target residue in a main-chain centered coord system, 
	 * appropriately rotated to match the original main-chain
	 * segment in the target. Average the structures in the target
	 * using the dodgy assumption that equivalent atoms are in the
	 * same order within the amino acid objects
	 */
	for (k=0; k<Ctmask.on_no(); k++)
	    Taa[k].set_values(0.0);
	for (hi=0; hi<Homs.len(); hi++)
	{
	    if (Equiv[hi]<0) continue;
	    
	    Aacid_& Haa=Homs[hi].Aas[Equiv[hi]];
	    for (k=0; k<Ctmask.on_no(); k++)
		Taa[k]+=Haa[k];	// may fail if ordering is silly, but...
	}
	Taa*=(1.0/Eqno);
	Taa+=Tctr;  // shift back to original position, with the new side chain
    }	    // for ti
}
// END of make_sidechain()

/* equiv_atoms(): switches off those atoms in Target which have
 * no active equivalents in Hom. 
 */
static void equiv_atoms(const Aacid_& Hom, Aacid_& Target)
{
    // same amino acid: all atoms are equivalent
    if (Hom.res_id()==Target.res_id())
    {
	Bits_ Tmask=Target.mask();
	Tmask&=Hom.mask();
	Target.mask(Tmask);
	return;
    }
    
    /* Make a local all-active target copy and switch all
     * atoms off: then switch the atoms equivalent to
     * the ones in Hom on, and finally (at the FINISH label)
     * re-mask the Target by AND-ing Eq's mask with its mask
     */
    Aacid_ Eq(Target);
    Eq.mask(false);

    // main chain atoms are always equivalent
    Eq.active("N", Hom.active("N"));
    Eq.active("CA", Hom.active("CA"));
    Eq.active("C", Hom.active("C"));
    Eq.active("O", Hom.active("O"));
    
    // Gly has no side chain
    if (Target.res_id()=='G' || Hom.res_id()=='G')
	goto FINISH;
    
    // C-beta: always equivalent (and that's all for Ala)
    Eq.active("CB", Hom.active("CB"));
    if (Target.res_id()=='A' || Hom.res_id()=='A')
	goto FINISH;
    
    // gamma: was there a beta-branch? (no equivalences among branched atoms)
    if (strchr("ITV", Hom.res_id()) || strchr("ITV", Target.res_id()))
	goto FINISH;
    
    // From now on, equivalent atom names may differ
    char Tname[4], Hname[4];
    
    // do the rest of the gamma position
    if (Target.res_id()=='C') strcpy(Tname, "SG");
    else if (Target.res_id()=='S') strcpy(Tname, "OG");
    else strcpy(Tname, "CG");
    if (Hom.res_id()=='C') strcpy(Hname, "SG");
    else if (Hom.res_id()=='S') strcpy(Hname, "OG");
    else strcpy(Hname, "CG");
    Eq.active(Tname, Hom.active(Hname));
    
    // end for those which cannot have equivalent delta positions
    if (strchr("CHPSW", Target.res_id()) || strchr("CHPSW", Hom.res_id()))
	goto FINISH;
    
    // aromatics: Phe,Tyr
    if (strchr("FY", Target.res_id()) && strchr("FY", Hom.res_id()))
    {
	Eq.active("CD1", Hom.active("CD1"));
	Eq.active("CD2", Hom.active("CD2"));
	Eq.active("CE1", Hom.active("CE1"));
	Eq.active("CE2", Hom.active("CE2"));
	Eq.active("CZ", Hom.active("CZ"));
    }
    if (strchr("FY", Target.res_id()) || strchr("FY", Hom.res_id()))
	goto FINISH;
    
    // delta-branch (cf. beta-branch above) D<->L matchable (both symmetric)
    if (Hom.res_id()=='L' && Target.res_id()=='D')
    {
	Eq.active("OD1", Hom.active("CD1"));
	Eq.active("OD2", Hom.active("CD2"));
	goto FINISH;
    }
    if (Hom.res_id()=='D' && Target.res_id()=='L')
    {
	Eq.active("CD1", Hom.active("OD1"));
	Eq.active("CD2", Hom.active("OD2"));
	goto FINISH;
    }
    if (strchr("DNL", Hom.res_id()) || strchr("DNL", Target.res_id()))
	goto FINISH;
    
    // unbranched delta remaining: EQKMR, all "CD" except M
    strcpy(Tname, (Target.res_id()=='M')? "SD": "CD");
    strcpy(Hname, (Hom.res_id()=='M')? "SD": "CD");
    Eq.active(Tname, Hom.active(Hname));
    
    // branched epsilon: don't map Q onto E
    if (strchr("EQ", Hom.res_id()) || strchr("EQ", Target.res_id()))
	goto FINISH;
    
    // rest of epsilon: KMR
    strcpy(Tname, (Target.res_id()=='R')? "NE": "CE");
    strcpy(Hname, (Hom.res_id()=='R')? "NE": "CE");
    Eq.active(Tname, Hom.active(Hname));
    if (Target.res_id()=='M' || Hom.res_id()=='M') goto FINISH;
    
    // zeta: K or R
    strcpy(Tname, (Target.res_id()=='R')? "CZ": "NZ");
    strcpy(Hname, (Hom.res_id()=='R')? "CZ": "NZ");
    Eq.active(Tname, Hom.active(Hname));
    
    // final cleanup: AND together the masks
    FINISH:
    Bits_ Tmask=Eq.mask();
    Tmask&=Target.mask();
    Target.mask(Tmask);
}
// END of equiv_atoms()

// ---- Output ----

/* pdb_list(): writes the sidechain-decorated target to a PDB
 * file Pdbfname. Only the coords and a few remarks are written.
 */
static void pdb_list(const char *Pdbfname, const char *Alignname, 
	const char *Modelname, const char *Homname, 
	const Array_<Aacid_>& Model)
{
    // create the entry, C style
    Pdbentry_ *Entry=(Pdbentry_ *) malloc(sizeof(Pdbentry_));
    
    // init some of the easier fields, zero others
    strcpy(Entry->Header, "PROTEIN MODEL");
    
    time_t Now=time(NULL);  // get current date
    strftime(Entry->Date, 10, "%d-%b-%y", localtime(&Now));
    
    strcpy(Entry->Pdbcode, "0XXX");
    Entry->Compound=(char *) calloc(61, sizeof(char));
    strcpy(Entry->Compound, "POLYPEPTIDE CHAIN");
    Entry->Source=(char *) calloc(61, sizeof(char));
    strcpy(Entry->Source, "DRAGON\'S SIDECHAIN HOMOLOGY MODELLER");
    strcpy(Entry->Expdta, "THEORETICAL MODEL");
    Entry->Resol=-1.0;
    Entry->Hbonds=NULL; Entry->Hbno=0;
    Entry->Ssbs=NULL; Entry->Ssbno=0;
    
    // create the chain
    Entry->Chainno=1;
    Entry->Chains=(Chain_*) malloc(sizeof(Chain_));
    
    // no secondary structure, H-bond and S-S bond list
    Entry->Chains->Secs=NULL; Entry->Chains->Secsno=0;
    Entry->Chains->Hbonds=NULL; Entry->Chains->Hbno=0;
    Entry->Chains->Ssbs=NULL; Entry->Chains->Ssbno=0;
    
    // generate the amino acid sequence
    int i, k, Rno=Model.len();
    Entry->Chains->Aano=Rno;
    Entry->Chains->Seq=(char *) calloc(Rno+1, sizeof(char));
    for (i=0; i<Rno; i++)
	Entry->Chains->Seq[i]=Model[i].res_id();
    Entry->Chains->Seq[Rno]='\0';   // para...
    
    // set chain ID and chain type
    Entry->Chains->Chid=' ';
    Entry->Chains->Type='P';	// "protein"
    
    // determine no. of atoms, alloc atom array
    int Atomno=0;
    for (i=0; i<Model.len(); i++)
	Atomno+=Model[i].active_len();
    Atom_ *A=(Atom_*) calloc(Atomno, sizeof(Atom_));
    
    // transfer atom names and coords
    const char *Atname=NULL;
    int p;
    const Vector_ *Cptr=NULL;
    for (i=k=0; i<Model.len(); i++)
    {
	for (p=0; p<Model[i].len(); p++)
	{
	    Atname=Model[i].name(p);
	    if (!Model[i].active(Atname)) continue;	// inactive
	    Cptr=Model[i].atom(Atname);
	    A[k].Atno=k+1;
	    strcpy(A[k].Id+1, Atname); A[k].Id[0]=' ';	// right justified
	    A[k].Alt=A[k].Rid=' ';
	    A[k].Aa=Entry->Chains->Seq[i];
	    A[k].Resno=i+1;
	    A[k].X=Cptr->operator[](0);
	    A[k].Y=Cptr->operator[](1);
	    A[k].Z=Cptr->operator[](2);
	    A[k].Occu=1.0; A[k].Bfact=0.0;
	    ++k;
	}
    }
    Entry->Chains->Atomno=Atomno;
    Entry->Chains->Atoms=A;
    
    // generate remarks
    const unsigned int REMARK_NO=3, REMARK_LEN=61;
    unsigned int ri;
    char **Remarks=new char* [REMARK_NO];
    for (ri=0; ri<REMARK_NO; ri++)
	Remarks[ri]=new char [REMARK_LEN];
    
    sprintf(Remarks[0], "ALIGNMENT FILE: %s", Alignname);
    sprintf(Remarks[1], "TARGET FILE: %s", Modelname);
    sprintf(Remarks[2], "HOMOLOGOUS STRUCTURES: %s", Homname);
    
    // write to disk (this is a C utility from "pdbprot") 
    put_pdb(Pdbfname, Entry, Remarks, REMARK_NO);
    
    // clean up, using the C utility fn from "pdbprot" for the malloc/free symmetry
    for (ri=0; ri<REMARK_NO; ri++)
	delete [] Remarks[ri];
    delete [] Remarks;
    free_pdb(Entry);
}
// END of pdb_list()

// ==== END OF PROGRAM sidech.c++ ====
