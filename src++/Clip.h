#ifndef CLIP_CLASS
#define CLIP_CLASS

// ==== PROJECT DRAGON: HEADER Clip.h ====

/* The command-line interpreter. */

// SGI C++ 6.2 (o32), IRIX 6.2, 21. Aug. 1996. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <fstream.h>

// ---- MODULE HEADERS ----

#include "Params.h"
#include "String.h"

// ---- TYPEDEFS ----

/* Runfnc_: a function ptr type that points to a function
 * that executes a run for a number of times 
 * and returns 0 if OK or a signal value
 * caught within. The function is supposed to catch SIGINT (Ctrl-C)
 * for interrupted runs.
 * Currently the function to be run by Clip_ is the dragon_run()
 * simulation routine (see Dragon.c++).
 */
typedef unsigned int (*Runfnc_)(unsigned int);

// ==== CLASSES ====

/* Class Clip_: a rudimentary command-line interpreter for DRAGON
 * which accepts commands from stdin (the default) or a script file.
 * All commands start with lowercase letters to distinguish them
 * from parameters (start with uppercase) which will be passed
 * to a Params_ object for interpretation.
 */
class Clip_
{
    // data 
    private:
    
    Params_& Params;	// a ref to the parameters
    String_ Prompt;	// the prompt in interactive mode
    int Cmdlevel;  // recursion depth (when a cmd file calls another)
    
    // methods
    public:
    
	// constructors
    /* Inits Clip so that it will be associated with
     * the parameter object Pars (not modifiable).
     * The prompt string is in Pr.
     */
    Clip_(Params_& Pars, const String_& Pr): Params(Pars), Prompt(Pr), Cmdlevel(0) {}
    
    /* get_command(): tries to obtain the next command from the file
     * named Cmd (cin if Cmd==NULL). Most commands
     * can be dealt with internally. If a "r[un] x"
     * command is seen which means "run run_func() x times
     * with the current parameters" then x>0 is returned. 
     * run_func() takes an unsigned int specifying the number of cycles to
     * be executed. It is supposed to return 0 at the end of a
     * successful run or the value of a signal caught inside.
     * If "q[uit]" is seen then -1 is returned, meaning "finish up everything".
     * If a script was processed and the end has been reached then
     * get_command() returns with 0. 
     */
    int get_command(const char *Cmd, Runfnc_ run_func);
    
    // private methods
    private:
    void put_prompt() const;

    // "forbidden methods"
    private:
    Clip_(const Clip_&);
    Clip_& operator=(const Clip_&);
    
};
// END OF CLASS Clip_

// ==== END OF HEADER Clip.h ====
#endif	/* CLIP_CLASS */
