// ==== PROJECT DRAGON: METHODS Params.c++ ====

/* Reads and writes DRAGON parameter files and keeps track of
 * the global parameters in between.
 */

// SGI C++ 7.1, IRIX 6.2, 15-Oct-1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <fstream.h>
#include <strstream.h>
#include <limits.h>
#include <float.h>

// ---- MODULE HEADER ----

#include "Params.h"

// ==== Params_ METHODS ====

// ---- Constructor ----

/* Inits all parameters to their default values. */
Params_::Params_():
	Strs(10), Longs(Paramlim_<long>(0, 0, LONG_MAX), 7), 
	Dbls(Paramlim_<double>(0.0, 0.0, DBL_MAX), 6)
{
    // parameter strings
    
    Strs[0].set_default("$DRAGON_DATA/DEFAULT.aln", 256);
    Strs[0].name_descr("Alnfnm", "Alignment file");
    
    Strs[1].set_default("$DRAGON_DATA/DEFAULT.pho", 256);
    Strs[1].name_descr("Phobfnm", "Amino acid hydrophobicity file");
    
    Strs[2].set_default("$DRAGON_DATA/DEFAULT.vol", 256);
    Strs[2].name_descr("Volfnm", "Side chain volume file");
    
    Strs[3].set_default("$DRAGON_DATA/DEFAULT.acd", 256);
    Strs[3].name_descr("Adistfnm", "File holding atom distances from C-alpha and sidechain centroids");
    
    Strs[4].set_default("$DRAGON_DATA/DEFAULT.sim", 256);
    Strs[4].name_descr("Simfnm", "Amino acid similarity matrix file");
    
    Strs[5].set_default("", 256);
    Strs[5].name_descr("Restrfnm", "External restraint file");
    
    Strs[6].set_default("", 256);
    Strs[6].name_descr("Sstrfnm", "Secondary structure assignment file");
    
    Strs[7].set_default("", 256);
    Strs[7].name_descr("Accfnm", "Surface/buried residue assignment file");
    
    Strs[8].set_default("", 256);
    Strs[8].name_descr("Homfnm", "Homologous structure PDB file");
    
    Strs[9].set_default("DRAGON_OUT", 256); // ".pdb" may be added in main
    Strs[9].name_descr("Outfnm", "Result PDB file");
    
    // integer parameters
    
    Longs[0].name_descr("Masterno", "Master sequence number (0=consensus)");

    Longs[1].set_deflims(40, 1, 500);
    Longs[1].name_descr("Maxiter", "Maximal number of iterations in 3D");
    
    Longs[2].name_descr("Randseed", "RNG seed");
    
    Longs[3].set_deflims(5, 1, 100);
    Longs[3].name_descr("Tangiter", "Maximal number of detangling iterations");
    
    Longs[4].set_deflims(0, 0, 1);
    Longs[4].name_descr("Graph", "Graphics off/on (SGI version only)");
    
    Longs[5].set_deflims(2, 2, LONG_MAX);
    Longs[5].name_descr("Minsepar", "Minimal sequential separation for homology restraints");
    
    Longs[6].set_deflims(30, 10, 100);
    Longs[6].name_descr("Speciter", "Maximal number of Specgrad optimisation iterations");
    
    // floating-point parameters
    
    Dbls[0].name_descr("Minscore", "Minimal score limit");
    Dbls[1].name_descr("Minchange", "Minimal relative score change");
    
    Dbls[2].set_deflims(0.999, 0.0, 1.0);
    Dbls[2].name_descr("Evfract", "Fraction of eigenvalues kept");
    
    Dbls[3].set_deflims(0.00636, 0.001, 0.012);
    Dbls[3].name_descr("Density", "Residue density [1/A^3]");
    
    Dbls[4].set_deflims(5.0, 0.0, DBL_MAX);
    Dbls[4].name_descr("Maxdist", "Maximal length of homology distance restraints");

    Dbls[5].set_deflims(0.02, 0.0001, 0.1);
    Dbls[5].name_descr("Speceps", "Precision for Specgrad iterations");
}

// ---- Access ----

/* reset_default(): resets all parameters to their default values. */
void Params_::reset_default()
{
    register unsigned int i;
    for (i=0; i<Strs.len(); i++)
	Strs[i].reset_default();
    
    for (i=0; i<Longs.len(); i++)
	Longs[i].reset_default();
    
    for (i=0; i<Dbls.len(); i++)
	Dbls[i].reset_default();
}
// END of reset_default()

/* changed(): returns the Changed value of the parameter called Parname.
 * These are set to true immediately after a successful input or default reset
 * and set to false after the first value read-out (see the x_value() methods
 * family below). The purpose of all this fuss is to save reconstruction
 * of big objects which may depend on global parameter values.
 */
bool Params_::changed(const String_& Parname) const
{
    register unsigned int i;
    for (i=0; i<Strs.len(); i++)
	if (Parname==Strs[i].name())	// found
	    return(Strs[i].changed());
    
    for (i=0; i<Longs.len(); i++)
	if (Parname==Longs[i].name())	// found
	    return(Longs[i].changed());
    
    for (i=0; i<Dbls.len(); i++)
	if (Parname==Dbls[i].name())	// found
	    return(Dbls[i].changed());
    
    // problem
    cerr<<"\n? Params_::changed("<<Parname<<") not found\n";
    return(false);
}
// END of changed()

/* reset_changed(): set the Changed bit of the parameter Parname
 * to false or of all parameters if Parname=="" (the default).
 * Return value: the number of bits flicked.
 */
int Params_::reset_changed(const String_& Parname)
{
    int Flick=0, All=(Parname==String_(""));
    
    register unsigned int i;
    for (i=0; i<Strs.len(); i++)
	if (All || Parname==Strs[i].name())	// found
	{
	    if (Strs[i].changed())
	    {
		Strs[i].not_changed();	// reset
		Flick++;
	    }
	    if (!All) return(Flick);
	}
    
    for (i=0; i<Longs.len(); i++)
	if (All || Parname==Longs[i].name())	// found
	{
	    if (Longs[i].changed())
	    {
		Longs[i].not_changed();	// reset
		Flick++;
	    }
	    if (!All) return(Flick);
	}
    
    for (i=0; i<Dbls.len(); i++)
	if (All || Parname==Dbls[i].name())	// found
	{
	    if (Dbls[i].changed())
	    {
		Dbls[i].not_changed();	// reset
		Flick++;
	    }
	    if (!All) return(Flick);
	}
    
    // problem
    if (!All)
	cerr<<"\n? Params_::reset_changed("<<Parname<<") not found\n";
    return(Flick);
}
// END of reset_changed()

/* s_value(), i_value(), f_value(): return the value of the parameter
 * called Parname, or NULL, 0, 0.0 if there was no such name in the
 * calling object (plus a warning is printed). These methods implement
 * a crude associative array but the type check is up to the programmer.
 */
const char* Params_::s_value(const String_& Parname)
{
    for (register unsigned int i=0; i<Strs.len(); i++)
	if (Parname==Strs[i].name())	// found
	{
	    Strs[i].not_changed();  // reset change status
	    return(Strs[i]);
	}
    
    // problem
    cerr<<"\n? Params_::s_value("<<Parname<<") not found, NULL returned\n";
    return(NULL);
}
// END of s_value()

long Params_::i_value(const String_& Parname)
{
    for (register unsigned int i=0; i<Longs.len(); i++)
	if (Parname==Longs[i].name())	// found
	{
	    Longs[i].not_changed();	// reset change status
	    return(Longs[i]);
	}
    
    // problem
    cerr<<"\n? Params_::i_value("<<Parname<<") not found, 0 returned\n";
    return(0);
}
// END of i_value()

double Params_::f_value(const String_& Parname)
{
    for (register unsigned int i=0; i<Dbls.len(); i++)
	if (Parname==Dbls[i].name())	// found
	{
	    Dbls[i].not_changed();  // reset change status
	    return(Dbls[i]);
	}
    
    // problem
    cerr<<"\n? Params_::f_value("<<Parname<<") not found, 0.0 returned\n";
    return(0.0);
}
// END of f_value()

// ---- Input/output ----

/* read_file(): reads parameter values from a file Fname.
 * Returns 1 on success, 0 on error.
 */
int Params_::read_file(const char *Fname)
{
    ifstream Inf(Fname);
    if (!Inf)
    {
	cerr<<"\n? Params_::read_file(\""<<Fname<<"\"): Cannot open\n";
	return(0);
    }
    Inf>>(*this);   // input
    Inf.close();
    return(1);	// OK
}
// END of read_file()

/* >>: reads parameter descriptions from a stream. The stream is
 * read line-by-line, with the following syntax:-
 * "NAME value \n".
 * Empty lines and lines beginning with '#' are considered comments
 * and skipped. Non-comment lines are read into a buffer and passed
 * to the read_from() methods of all parameter objects contained
 * by the Params_ object being updated. This way, only the parameter
 * with the correct name will be updated.
 */
istream& operator>>(istream& In, Params_& P)
{
    const unsigned int LINELEN=132;
    char Line[LINELEN];
    istrstream Instr(Line, LINELEN);

    while (In.good())
    {
        memset(Line, '\0', LINELEN);
	In.getline(Line, LINELEN, '\n');
	if (!strlen(Line) || Line[0]=='#')  // skip empty and comment
	    continue;
	    
	/* Pass the input buffer to all member objects of P.
	 * The read_from() method returns 0 if the line was not
	 * consumed and therefore can be passed on to another
	 * object. 1 is returned if OK, -1 if an input error
	 * occurred: both non-0 conditions mean that no further
	 * processing is necessary and we can 'continue'.
	 */
	Instr.seekg(ios::beg); Instr.clear();
	unsigned int i;
	int Rdstat=0;
	
	// read strings first
	for (i=0; !Rdstat && i<P.Strs.len(); i++)
	{
	    Rdstat=P.Strs[i].read_from(Instr);
	}
	if (Rdstat) continue;

	// read integral values
	for (i=0; !Rdstat && i<P.Longs.len(); i++)
	{
	    Rdstat=P.Longs[i].read_from(Instr);
	}
	if (Rdstat) continue;

	// read floating-point
	for (i=0; !Rdstat && i<P.Dbls.len(); i++)
	{
	    Rdstat=P.Dbls[i].read_from(Instr);
	}
    }
    return(In);
}
// END of >>

/* write_file(): lists the complete parameter list to a file Fname.
 * Returns 0 on error, 1 if OK.
 */
int Params_::write_file(const char *Fname) const
{
    ofstream Outf(Fname);
    if (!Outf)
    {
	cerr<<"\n? Params_::write_file(\""<<Fname<<"\"): Cannot open\n";
	return(0);
    }
    Outf<<(*this);   // output
    Outf.close();
    return(1);	// OK
}
// END of write_file()

/* <<: prints parameter descriptions to a stream Out. */
ostream& operator<<(ostream& Out, const Params_& P)
{
    unsigned int i;
    
    Out<<"# --- String parameters ---\n";
    for (i=0; i<P.Strs.len(); i++)
	P.Strs[i].write_to(Out);
        
    Out<<"\n# --- Integer-valued parameters ----\n";
    for (i=0; i<P.Longs.len(); i++)
	P.Longs[i].write_to(Out);
        
    Out<<"\n# ---- Floating-point parameters ----\n";
    for (i=0; i<P.Dbls.len(); i++)
	P.Dbls[i].write_to(Out);
        
    return(Out);
}
// END of <<

/* list_changed(): lists all parameters that have been changed
 * to the output stream Out (default cout). No comments are printed.
 * Return value: the number of changed parameters.
 */
unsigned int Params_::list_changed(ostream& Out) const
{
    unsigned int i, Chgno=0;

    for (i=0; i<Strs.len(); i++)
	if (Strs[i].changed())
	    { Strs[i].write_to(Out, false); Chgno++; }
        
    for (i=0; i<Longs.len(); i++)
	if (Longs[i].changed())
	    { Longs[i].write_to(Out, false); Chgno++; }
        
    for (i=0; i<Dbls.len(); i++)
	if (Dbls[i].changed())
	    { Dbls[i].write_to(Out, false); Chgno++; }
        
    return(Chgno);
}
// END of list_changed()

/* list_param(): lists the parameter called Parname to the output stream
 * Out (default cout). Prints a message to cerr and returns 0 if 
 * there was no such parameter, otherwise returns 1.
 */
int Params_::list_param(const String_& Parname, ostream& Out) const
{
    register unsigned int i;
    for (i=0; i<Strs.len(); i++)
	if (Parname==Strs[i].name())	// found
	{
	    Strs[i].write_to(Out);
	    return(1);
	}
    
    for (i=0; i<Longs.len(); i++)
	if (Parname==Longs[i].name())	// found
	{
	    Longs[i].write_to(Out);
	    return(1);
	}
    
    for (i=0; i<Dbls.len(); i++)
	if (Parname==Dbls[i].name())	// found
	{
	    Dbls[i].write_to(Out);
	    return(1);
	}
    
    // problem
    cerr<<"\n? Params_::list_param("<<Parname<<") not found\n";
    return(0);
}
// END of list_param()

// ==== END OF METHODS Params.c++ ====
