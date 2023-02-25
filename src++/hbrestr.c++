// ==== PROGRAM hbrestr.c++ ====

/* Generates H-bond distance restraints for secondary
 * structures in QUANTA format.
 */

// SGI C++, IRIX 5.3, 26. Apr. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>

// ---- MODULE AND UTILITY HEADERS ----

#include "Pieces.h"	// secstr storage
#include "pdbprot.h"	// PDB I/O and storage (in C)
#include "Vector.h"	// vector algebra

// ---- PROTOTYPES ----

static void write_restraints(const Helix_& Hel, ostream& Out);

static void print_restraint(int Co, int Nh, ostream& Out);

// ==== MAIN ====

/* The program takes three command-line arguments: a PDB filename
 * and a DRAGON-IV secondary structure description file name
 * and an output file name.
 * The secstr-s must fit into the range of residues in the PDB file
 * but otherwise no checks are made to ensure that the assignments
 * correspond to the actual structure.
 * Output is a QUANTA (4.x) *.dcn "NOE" distance restraint file format
 * to the third argument. 
 */
int main(int argc, char *argv[])
{
    // three arguments must be present
    if (argc<4)
    {
	cerr<<"\n! Usage: "<<argv[0]<<" PDB_file secstr_file constraint_file\n";
	exit(EXIT_FAILURE);
    }
    
    Pdbentry_ *Pdb=get_pdb(argv[1], ALLATOMS, RELAXED);
    if (Pdb==NULL || !Pdb->Chainno)
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read PDB file \""<<argv[1]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    
    // use the first chain only
    Chain_ *Chain=Pdb->Chains;
    int Resno=Chain->Aano;
    Pieces_ Pieces(Resno);  // init to hold Resno residues
    
    // get the secstr assignment
    if (!Pieces.read_secstr(argv[2]))
    {
	cerr<<"\n! "<<argv[0]<<": Cannot read secstr file \""<<argv[2]<<"\"\n";
	exit(EXIT_FAILURE);
    }
    
    // open the output file and print the prologue
    ofstream Out(argv[3]);
    Out<<"*CHARMm distance constraints faked by "<<argv[0]<<endl;
    Out<<"*PDB file: \""<<argv[1]<<"\"\n";
    Out<<"*DRAGON-IV secstr file: \""<<argv[2]<<"\"\n";
    Out<<"NOE\nRESET\n";
    
    // process the secstr list
    Clist1_<Sstr_> Slist=Pieces.secs();
    for (Slist.begin(); Slist!=NULL; Slist++)
    {
	if ((*Slist)->is_helix())
	    write_restraints(Helix_(*Slist), Out);	// helical version
	else
	    write_restraints(Beta_(*Slist), Chain, Out);	// sheet version
    }

    // finish up
    Out<<"SCALE     1.0000\nEND\n";    // epilogue
    Out.close();
    exit(EXIT_SUCCESS);
}

// ==== FUNCTIONS ====

/* write_restraints(): lists the main-chain >C=O...H-N< bond
 * restraints to the stream Out. There are versions for
 * helices and sheets: the latter needs the corresponding coordinates
 * (Chain) to find out about the H-bond pattern (depends on
 * which C-betas are "up" or "down", rather ugly...).
 */
static void write_restraints(const Helix_& Hel, ostream& Out)
{
    // get the helix phasing
    int Phase=4;    // for alpha or default
    if (Hel.helix_type()==Helix_::HX310) Phase=3;   // 3/10 helix
    else if (Hel.helix_type()==Helix_::HXPI) Phase=5;	// pi-helix
    
    /* trundle along the helix: the >C=O of the i:th residue
     * is always bound to the H-N< of the (i+Phase):th residue
     * which simplifies matters to a large extent
     */
    for (int i=Hel.beg(); i<=Hel.end()-Phase; i++)
	print_restraint(i, i+Phase, Out);
}

static void write_restraints(const Sheet_& Sh, const Chain_ *Chain, 
	ostream& Out)
{
    // I HATE THIS!!!! I ALWAYS HATED BETA-SHEETS AND
    // I HATE THE BLOODY BOOKKEEPING ASSOCIATED WITH THEM!!!!!
    // WILL NOT WRITE THIS BLOODY PROGRAM!!!!!!!!!!!!!
}

/* print_restraint(): prints the H-bond restraint to Out.
 * Co and Nh are residue numbers. WARNING: the string "0XXX"
 * gets printed instead of the QUANTA segment identifier!
 */
static void print_restraint(int Co, int Nh, ostream& Out)
{
    Out<<"ASSIGN SELE ATOM 0XXX "<<Co<<" O\t\tEND -\n";
    Out<<"       SELE ATOM 0XXX "<<Nh<<" HN\t\tEND -\n";
    Out<<"  KMIN   25.00 RMIN   1.900 KMAX   25.00 RMAX    2.10\n";
}
// END of print_restraint()

// ==== END OF PROGRAM hbrestr.c++ ====
