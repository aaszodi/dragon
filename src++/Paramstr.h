#ifndef PARAMSTR_HEADER
#define PARAMSTR_HEADER

// ==== PROJECT DRAGON: HEADER Paramstr.h ====

/* Stores global string (filename) parameters. */

// SGI C++ 7.1, IRIX 6.2, 7. Mar. 1997. Andris Aszodi

// ---- UTILITY HEADERS ----

#include "Parambase.h"

// ==== CLASSES ====

/* Class Paramstr_: stores a string parameter. Can be queried to return
 * this parameter which can be set via the input routine only. 
 * Derived from Parambase_
 */
class Paramstr_: public Parambase_
{
    // data
    protected:
    String_ Value;  // the string parameter
    String_ Default;	// the default string
    
    // methods
    public:
    
	// constructors
    /* Sets the maximal length of the string variables within to Size.
     * If Size==0 (default), then the Defval string's length is
     * used, if Size is less than the maximal length of Defval,
     * then it is adjusted.
     */
    Paramstr_(const char *Defval=NULL, unsigned int Size=0, 
		const char *Nm=NULL, const char *Ds=NULL);
    
	// access
    /* Converts the calling object to a (const) char ptr.
     * If the value of the calling object begins with the string
     * "$DRAGON_DATA/", then this string will be substituted
     * with the value of the environment variable DRAGON_DATA
     * if defined, otherwise it is removed.
     */
    operator const char* () const;
    
    /* resets the value to its default */
    void reset_default() { Value=Default; Changed=true; }
    
    /* set_default(): sets the default string to Defval and the maximal
     * size to Size (default==0). If Size==0, then the length of Defval
     * is used as the maximal size (same as in the ctor). Sets the
     * Value string to Defval also.
     * Returns maximal size.
     */
    unsigned int set_default(const char *Defval, unsigned int Size=0);
    
	// input/output
    
    /* read_from(): reads the following line from an input stream In:
     * "NAME string"
     * where NAME should match the Name member of the object and "string"
     * should be a string not containing whitespaces. 
     * If the line was of the form "NAME " then the default string is
     * used as the value (which would have been "" otherwise).
     * If "string" is longer than the max. size in the object, then
     * the extra chars are left in In but this is OK.
     * Return values: -1 on error, 0 if NAME did not match, 1 if OK.
     */
    int read_from(istream& In);
    
    /* Lists to Out. */
    ostream& write_to(ostream& Out, bool Comments=true) const;

};
// END OF CLASS Paramstr_

// ==== END OF HEADER Paramstr.h ====

#endif	/* PARAMSTR_HEADER */
