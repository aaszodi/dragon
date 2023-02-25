#ifndef ACDIST_CLASSES
#define ACDIST_CLASSES

// ==== PROJECT DRAGON: HEADER Acdist.h ====

/* For the storage and retrieval of side-chain atom
 * distances from the C-alpha atom or the sidechain centroid.
 */

// SGI C++, IRIX 6.2, 8. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <ctype.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "String.h"
#include "Array.h"

// ==== CLASSES ====

/* Class Acs_: this class stores the side-chain atom
 * distances from the C-alpha atom and the sidechain centroid
 * for each "standard" amino acid. Can be queried for these
 * distances. 
 */
class Acs_
{
    friend class Acdist_;
    
    /* struct Acd_: stores the name of the sidechain atom
     * according to the PDB conventions, and the distance of
     * the atom from the C-alpha and the sidechain centroid.
     */
    struct Acd_
    {
	String_ Atname;	    // atom name
	float Adist, Cdist; // CA and centroid distances
	
	Acd_(const char *Atn=NULL, float Ad=0.0, float Cd=0.0)
	    : Atname(Atn), Adist(Ad), Cdist(Cd) {}
    };
    // END OF STRUCT Acd_
    
    // data
    
    char Aa;	// the amino acid type
    Array_<Acd_> Acds;  // array of Acd_ entries
    
    // all methods are private
    
	// access
    /* ca_dist(), scc_dist(): return the distance of Atom from
     * the C-alpha or the sidechain centroid, respectively.
     * -1.0 returned and a warning printed if Atom was not found.
     */
    float ca_dist(const String_& Atom) const;
    float scc_dist(const String_& Atom) const;
    
    /* set_acd(): stores the distances of Atom from the CA and the
     * side chain centroid (Ad, Cd). Prints a warning and does
     * nothing if Atom was invalid or if Ad or Cd were negative.
     * Return value: 0 on error, non-0 otherwise.
     */
    int set_acd(const String_& Atom, float Ad, float Cd);
    
    /* get_acdptr(): returns a const pointer to the Acd_ record
     * for the side chain atom Atom if it is part of the side chain
     * of the calling object; otherwise, NULL is returned and a
     * warning is printed.
     */
    const Acd_* get_acdptr(const String_& Atom) const
    {
        if (!Acds.len()) // 12-Jul-1997
	  {
	    cerr<<"\n? Acs_::get_acdptr(): Empty amino acid atom list\n";
	    return(NULL);
	  }
	for (register unsigned int i=0; i<Acds.len(); i++)
	    if (Acds[i].Atname==Atom) return(&(Acds[i]));
	cerr<<"\n? Acs_::get_acdptr(): Atom \""<<Atom<<"\" not in side chain\n";
	return(NULL);
    }

	// setup
    /* set_dists(): sets up the sidechain atom distances (from CA and 
     * centroid) for amino acid Aac. If Aac is lowercase, then it
     * is converted to uppercase. If Aac is illegal, then it is
     * assumed to be 'X' ("anything" or "unknown"). 'X' will store
     * nothing at all, though. In all valid cases, an internal array
     * is set up to hold the names of the sidechain atoms (according to
     * the PDB convention) and their distances. By default, these
     * distances will be taken from the most abundant rotamers (>=10%) in
     * the Ponder/Richards library as defined in Quanta 4.1.
     */
    void set_dists(char Aac);
    
    void ala();
    void cys();
    void asp();
    void glu();
    void phe();
    void gly();
    void his();
    void ile();
    void lys();
    void leu();
    void met();
    void asn();
    void pro();
    void gln();
    void arg();
    void ser();
    void thr();
    void val();
    void trp();
    void tyr();
    void unk();
};
// END OF CLASS Acs_

/* Class Acdist_: this class stores an array of Acs_ objects,
 * one for each amino acid. The side-chain atom distances
 * from the C-alpha and sidechain centroid can be queried
 * safely here. Distance data may be read from a file.
 * This class represents the public interface to this module.
 */
class Acdist_
{
    // data
    private:
    
    static const char *AAcodes;   // legal amino acid codes (1-letter)
    Acs_ Acss[20];	    // simple array holding Acs_ objects for each amino acid
    
    // methods
    public:
    
	// constructor
    /* Inits to the distances in file Fname or to the defaults
     * (cf. Acs_::set_dists()) if Fname is "" or NULL (the default).
     */
    Acdist_(const char *Fname=NULL)
    {
	reset();
	read_file(Fname);
    }
    
	// access
    /* reset(): resets all amino acid sidechain distances to their default values
     * as defined in Acs_::set_dists(). 
     */
    void reset()
    {
	for (register unsigned int i=0; i<strlen(AAcodes); i++)
	  {
	    Acss[i].set_dists(AAcodes[i]);
	  }
    }
    
    /* ca_dist(), scc_dist(): return the distance of Atom from
     * the C-alpha or the sidechain centroid for amino acid Aa, respectively.
     * -1.0 returned and a warning printed if Aa or Atom was not found.
     */
    float ca_dist(char Aa, const String_& Atom) const;
    float scc_dist(char Aa, const String_& Atom) const;
    
	// input
    /* read_file(): reads sidechain atom distance data from 
     * the file called Fname. Does nothing if Fname is "", NULL (the default)
     * or cannot be opened. Updates only the distances which
     * are explicitly mentioned in the file; for a complete reset, 
     * call reset(). The format of a line in the distance file
     * is:-
     * "AAcode Atomname CAdist CTRdist\n"
     * where AAcode is a 1-letter amino acid code (char), Atomname is
     * an all-uppercase string (must be a PDB-type sidechain atom), 
     * CAdist is the distance of the atom from the C-alpha atom in
     * angstroms (float),  CTRdist is the distance of the atom from
     * the sidechain centroid (float). The items are separated by
     * whitespaces. Lines beginning with '#' are comments. Invalid
     * data elicit warnings and will be skipped.
     * Return value: 0 on error, 1 if OK.
     */
    int read_file(const char *Fname);
    
    /* >>: reads from the stream Inf. See comments to read_file() */
    friend istream& operator>>(istream& Inf, Acdist_& Acdist);
    
    // private and forbidden methods
    private:
    
    static int get_idx(char Aa)
    {
        if (!Aa || !isalpha(Aa)) return(-1);  // \0 was found by strchr() GCC/Linux
	Aa=toupper(Aa);
	char *Pos=strchr(AAcodes, Aa);
	return((Pos==NULL)? -1: Pos-AAcodes);
    }
    
    int set_acd(unsigned int Idx, const String_& Atom, float Ad, float Cd)
    { return(Acss[Idx].set_acd(Atom, Ad, Cd)); }

    Acdist_(const Acdist_&);	// no copy or assignment
    Acdist_& operator=(const Acdist_&);
};
// END OF CLASS Acdist_

// ==== END OF HEADER Acdist.h ====
#endif	/* ACDIST_CLASSES */

