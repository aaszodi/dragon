#ifndef PARAMS_HEADER
#define PARAMS_HEADER

// ==== PROJECT DRAGON: HEADER Params.h ====

/* Reads and writes DRAGON parameter files and keeps track of
 * the global parameters in between.
 */

// SGI C++ 7.1, IRIX 6.2, 7. Mar. 1997. Andris Aszodi

// ---- MODULE HEADERS ----

#include "Paramstr.h"
#include "Paramlim.h"

// ---- UTILITY HEADERS ----

#include "Array.h"

// ==== CLASSES ====

/* Class Params_: stores all parameters for DRAGON. Can read/write
 * a full parameter file or just a few lines from a stream.
 * Can be queried for the values of the parameters.
 */
class Params_
{
    // data
    private:
    
    Array_<Paramstr_> Strs;	// string parameters
    Array_<Paramlim_<long> > Longs;		// integer parameters
    Array_<Paramlim_<double> > Dbls;	// floating-point parameters
    
    // methods
    public:
   
	// constructors
    /* Inits all parameters to their default values. */
    Params_();
     
	// access

    /* reset_default(): resets all parameters to their default values. */
    void reset_default();
    
    /* changed(): returns the Changed value of the parameter called Parname.
     * These are set to true immediately after a successful input or default reset
     * and set to false after the first value read-out (see the x_value() methods
     * family below). The purpose of all this fuss is to save reconstruction
     * of big objects which may depend on global parameter values.
     */
    bool changed(const String_& Parname) const;
    
    /* reset_changed(): set the Changed bit of the parameter Parname
     * to false or of all parameters if Parname=="" (the default).
     * Return value: the number of bits flicked.
     */
    int reset_changed(const String_& Parname="");
    
    /* s_value(), i_value(), f_value(): return the value of the parameter
     * called Parname, or NULL, 0, 0.0 if there was no such name in the
     * calling object (plus a warning is printed). These methods implement
     * a crude associative array but the type check is up to the programmer.
     */
    const char* s_value(const String_& Parname);
    long i_value(const String_& Parname);
    double f_value(const String_& Parname);

	// input/output 
    
    /* read_file(): reads parameter values from a file Fname.
     * Returns 1 on success, 0 on error.
     */
    int read_file(const char *Fname);
    
    /* >>: reads parameter descriptions from a stream. The stream is
     * read line-by-line, with the following syntax:-
     * "NAME value \n".
     * Empty lines and lines beginning with '#' are considered comments
     * and skipped. Non-comment lines are read into a buffer and passed
     * to the read_from() methods of all parameter objects contained
     * by the Params_ object being updated. This way, only the parameter
     * with the correct name will be updated.
     */
    friend istream& operator>>(istream& In, Params_& P);
    
    /* write_file(): lists the complete parameter list to a file Fname.
     * Returns 0 on error, 1 if OK.
     */
    int write_file(const char *Fname) const;
    
    /* <<: prints parameter descriptions to a stream Out. */
    friend ostream& operator<<(ostream& Out, const Params_& P);

    /* list_changed(): lists all parameters that have been changed
     * to the output stream Out (default cout). No comments are printed.
     * Return value: the number of changed parameters.
     */
    unsigned int list_changed(ostream& Out=cout) const;
    
    /* list_param(): lists the parameter called Parname to the output stream
     * Out (default cout). Prints a message to cerr and returns 0 if 
     * there was no such parameter, otherwise returns 1.
     */
    int list_param(const String_& Parname, ostream& Out=cout) const;
};
// END OF CLASS Params_

// ==== END OF HEADER Params.h ====

#endif	/* PARAMS_HEADER */
