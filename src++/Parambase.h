#ifndef PARAMBASE_HEADER
#define PARAMBASE_HEADER

// ==== PROJECT DRAGON: HEADER Parambase.h ====

/* Abstract base class for global parameter storage.
 * There is no associated Parambase.c++ file.
 */

// SGI C++ 7.1, IRIX 6.2, 7. Mar. 1997. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "String.h"

/* NOTE: The N32/N64 compilers recognise the built-in bool type
 * and define the macro _BOOL. This workaround is provided for
 * the O32 compiler.
 */
#if !(defined(_BOOL) || defined(HAS_SEEN_BOOL))
#define HAS_SEEN_BOOL
typedef unsigned char bool;
static const bool false=0;
static const bool true=1;
#endif	/* _BOOL */

// ==== CLASSES ====

/* NOTE: Parameters are stored in a little class hierarchy:-
 * 
 *       [ Parambase_ ]
 *           |
 *     +-----+-------+
 *     |             |
 *     V             V
 * Paramstr_    Paramlim_<template>
 * 
 * where Parambase_ is an ABC. String parameters are in Paramstr_, 
 * the numerical parameters are in Paramlim_ template classes.
 */
 
/* Class Parambase_: abstract base class for the parameter storage.
 * Stores a "Name" string (the name of the parameter) and a
 * "Description" string (which explains briefly what the parameter is
 * and is printed as a comment).
 */
class Parambase_
{
    // data
    protected:
    
    String_ Name;   // the parameter name
    String_ Descr;  // the description
    bool Changed;    // false if reset, true after input
    
    // methods
    public:
    
	// constructors and destructor
    /* Inits Name to Nm and Descr to Ds (both NULL by default which
     * sets the String_ variables to "").
     */
    Parambase_(const char *Nm=NULL, const char *Ds=NULL):
	Name(Nm), Descr(Ds), Changed(true) {}
    virtual ~Parambase_() {}
    
	// access
    /* Resets the parameter to its default value. */
    virtual void reset_default() =0;
    
    /* changed(): returns Changed (true after input, false if no change) */
    bool changed() const { return(Changed); }
    
    /* not_changed(): resets Changed to false. Used after value enquiries. */
    void not_changed() { Changed=false; }
    
    /* name(): returns the name. */
    const String_& name() const { return(Name); }
    
    /* name_descr(Ds): sets the name and the description. */
    void name_descr(const char *Nm, const char *Ds) { Name=String_(Nm), Descr=String_(Ds); }
    
	// input/output
    
    /* read_from(): reads the following line from an input stream In:
     * "NAME value"
     * where NAME should match the Name member of the object and value
     * should conform to the actual type of Value. If NAME does not
     * match then Value will not be changed, the input stream will be reset to
     * its previous position and the method returns w/o an error. This
     * mechanism serves to skip other irrelevant parameter lines. If NAME
     * matches then input is attempted: if it is unsuccessful, then
     * the fail bit is set, Value will be Default and a warning is printed to cerr.
     * Misformed lines are consumed.
     * Return values: -1 on error, 0 if NAME did not match, 1 if OK.
     * Note that the non-0 return means that the line was consumed, 
     * 0 means that the input line was intended for another object and
     * is reset for further processing by the others (cf. >> in Params_).
     */
    virtual int read_from(istream& In) =0;
    
    /* write_to(): writes the calling object to the output stream Out
     * in the following format if Comments==true (the default):-
     * "# DESCRIPTION_STRING (default DEFAULT)\n"
     * "NAME value\n"
     * which can be parsed by the input routine. If Comments==false, 
     * then the "#" line is omitted.
     */
    virtual ostream& write_to(ostream& Out, bool Comments=true) const =0;
};
// END OF CLASS Parambase_

// ==== END OF HEADER Parambase.h ====

#endif	/* PARAMBASE_HEADER */
