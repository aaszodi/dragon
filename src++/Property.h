#ifndef PROPERTY_CLASS
#define PROPERTY_CLASS

// ==== PROJECT DRAGON: HEADER Monomer.h ====

/* Amino-acid property data storage. */

// SGI C++ 4.0, IRIX 5.3, 12. Jan. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>

// ---- DEFINITIONS ----

#define AANO 24	    /* no. of amino acids: 20 + B,Z,X */

// ---- TYPEDEFS ----

/* amino acid property array type */
typedef float Proparray_[AANO];

// ==== CLASSES ====

/* Class Property_: stores the list of the 1-letter amino acid codes
 * as a const static string and stores an AANO-long traditional C
 * float array of some amino acid property. 
 */
class Property_
{
    private:
    
    // data
    static const char GAP;
    static const char Aacodes[AANO+1];	// amino acid codes + GAP
    
    Proparray_ Prop;	// the property array
    const float *Default;   // the corresponding default
    
    public:

    // default properties
    static const Proparray_ Hyphobdef;	// Levitt's hydrophobicity
    static const Proparray_ Volumedef;	// Amino acid side chain volumes
    static const Proparray_ Cabdistdef;	// C-alpha:sidechain ctr dist (not squared)
    
    // methods
    
	// constructor
    /* Inits so that the default values are taken from the external
     * location pointed to by Defptr. Prop is initialised to { 0.0 ... }
     * if Defptr==NULL (the default ;-)). Defptr is assumed to
     * point to something meaningful such as the static const arrays
     * provided here.
     */
    Property_(const float *Defptr=NULL);
    
	// access
    /* reset(): resets to the default values in the array pointed to
     * by Default, or to {0.0 ...} if it was NULL (see ctor above).
     */
    void reset();
    
    /* [X]: returns the property value for amino acid Aa. If Aa is
     * lowercase, then it will be converted to uppercase. Invalid
     * chars will be treated as 'X' (unknown) silently.
     */
    float operator[](char Aa) const;
    
    /* avg_val(): returns the average value of the property for the
     * amino acids in the string Posstr (representing an alignment position).
     */
    float avg_val(const char *Posstr) const;
    
    /* read_file(): reads the property data file Fname. The format is:-
     * <char> <value> ...\n
     * where <char> must be a valid uppercase 1-letter amino acid code, 
     * and <value> is a floating-point number. Lines beginning with '#'
     * are comments and ignored. Non-parsable lines will be skipped
     * with a warning to cerr. 
     * NOTE: the format has changed from Version 3.x and the "title line"
     * is no longer supported. Put a '#' in front of the titles in the
     * old property data files to avoid error messages.
     * Return value: 0 on error, 1 if OK.
     */
    int read_file(const char *Fname);
    
    /* >>: read from stream In. See comments to read_file(). */
    friend istream& operator>>(istream& In, Property_& Pr);
    
    /* Output: list to stream Out neatly. */
    friend ostream& operator<<(ostream& Out, const Property_& Pr);
    
    // private methods
    private:
    int get_index(char Aa) const;
    
    // "forbidden methods": no copying or assignment
    Property_(const Property_&);
    Property_& operator=(const Property_&);
    
};
// END OF CLASS Property_

// ==== END OF HEADER Property.h ====

#endif	/* PROPERTY_CLASS */
