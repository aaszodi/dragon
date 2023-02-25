#ifndef SIMIL_CLASS
#define SIMIL_CLASS

// ==== PROJECT DRAGON: HEADER Simil.h ====

/* Implements a class for amino acid similarity matrices. */

// SGI C++ 4.0, IRIX 5.3, 12. Jan. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "Sqmat.h"

// ==== CLASSES ====

/* Class Simil_ : amino acid similarity class. Holds a square
 * matrix for similarity values and an amino acid permutation string
 * since most similarity matrices in the literature do not come in
 * the simple "ABC..." ordering. Can be queried for individual
 * amino acid (or gap) pair similarities and can calculate the
 * normalised similarity measure for a set of amino acids.
 * The matrices are shifted so that all entries are non-negative.
 * The internal values can be changed only via reading a suitably formatted
 * disk file; see the read_file() method.
 */
class Simil_
{
    // data
    private:
    
    char *Aacodes;  // 1-letter amino acid codes in a string
    Sqmat_ Sim;	    // the similarity matrix
    
    // methods
    public:
    
	// constructors
    /* Inits to the 26 uppercase letters of the alphabet and to the unit matrix. */
    Simil_();
    
	// destructor
    ~Simil_() { delete [] Aacodes; }
    
	// access
    /* simil(): returns the similarity value between amino acids A and B.
     * Lowercase letters will be silently converted to uppercase, unknown
     * characters trigger a warning and will be converted to 'X'. If X is
     * not in the code string then 0.0 is returned.
     */
    double simil(char A, char B) const;

    /* cons(): calculates the normalised consensus value for the amino acid
     * codes in the string Aas. Empty strings have a value of 0.0, strings
     * with just 1 char (possibly coming from 1-sequence alignments) give
     * 1.0, "AAA..A" have also 1.0, anybody else between 0 and 1. Gaps are
     * skipped in the summation (indicated by '-'). Returns the consensus
     * character: the most abundant amino acid in the string. The consensus
     * value is returned in Consval. (The consensus char is not always needed.)
     */
    char cons(const char *Aas, double& Consval) const;
    
    /* reset(): brings the object back to the same state it was in when
     * the default ctor created it (26 uppercase letters as code string, 
     * unit matrix as similarity matrix).
     */
    void reset();
    
	// input/output
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
    int read_file(const char *Fname);
    
    /* >>: reads a similarity matrix from Inf. See comments to read_file() */
    friend istream& operator>>(istream& Inf, Simil_& Simil);
    
    friend ostream& operator<<(ostream& Out, const Simil_& S);
    
    // private methods
    private:

    int pos(char& C) const;
    
    // "forbidden functions": no copy or assignment
    Simil_(const Simil_&);
    Simil_& operator=(const Simil_&);
    
};
// END OF CLASS Simil_

// ==== END OF HEADER Simil.h ====

#endif	/* SIMIL_CLASS */

