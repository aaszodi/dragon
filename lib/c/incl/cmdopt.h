#ifndef CMDOPT_H
#define CMDOPT_H

/* One-letter command line option processing. */

/* ANSI C, IRIX 5.2, 4. Apr. 1995. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* parse_optstr(): constructs and returns the hidden Option_ array from the
 * string Cmdoptstr. The length of the array is also stored in a hidden location.
 * The option string is composed of the following tokens separated
 * by whitespaces:-
 * 
 * "x": stands for the Boolean command line option -x
 * "xYz" : the command line options -x -Y -z (may be grouped together)
 * "x%d<name>" : option -x expecting a mandatory integer argument
 *               which is described as "name" in the help string.
 * "x%f<name>" : as above, expecting a floating-point argument
 * "x%s<name>" : as above, expecting a string argument
 * 
 * 'x' can be the characters 'a'..'z', 'A'..'Z', '0'..'9', '#'.
 * Invalid options are ignored. If nothing was found then Cmdopts==NULL and 
 * Cmdoptno==0. To be called only once.
 */
void parse_optstr(char *Cmdoptstr);

/* get_options(): processes the command line (argc,argv) to find options
 * encoded in the hidden array Cmdopts (length Cmdoptno). If an option is found then
 * the corresponding Cmdopts[] member is updated to show which argv[] member
 * contained it and if it had an argument, then it is attempted to be
 * stored as well. Options not found or having malformed arguments will
 * be set as "absent".
 * Return value: the argv[] index of the first non-option member.
 * The index is multiplied by -1 if an error occurred.
 */
int get_options(int argc, char *const *argv);

/* optval_[bool,int,dbl,str](): these functions query the hidden Cmdopts[] array
 * after the command line has been processed. Och is the option char
 * and the value is returned in *Val. If the option char is invalid
 * then a warning is printed. If the option was not present on the
 * command line then *Val is not updated and the functions return 0, 
 * otherwise the Idx field from the corresponding Cmdopts record is returned.
 * Val may be set to NULL.
 */
int optval_bool(char Och);
int optval_int(char Och, int *Val);
int optval_dbl(char Och, double *Val);
int optval_str(char Och, char **Val);

/* opt_defval(): this routine provides an alternative to the optval_[...]
 * functions. Och is the option character, Val is a void pointer to
 * the variable which should hold the option value, Defval is a
 * void pointer to the default value. If the option is found, 
 * then *Val is set to its value if Val!=NULL. If not found or on error, 
 * then *Val will be set to *Defval (to the default) if both
 * Val and Defval !=NULL. 
 * Be careful, there is NO POINTER TYPE CHECKING!
 * Return value: 0 on error or if the option was not found, 
 * the option array index otherwise.
 */
int opt_defval(char Och, void *Val, const void *Defval);

/* opt_helpstr(): generates a "help string" from the option list hidden in
 * Cmdopts[]. Boolean options are collected together, the "argumented"
 * options are listed separately. If -x and -y are switches, -i expects
 * an integer option and -D a double, then the following string will
 * be returned: "[-xy] [-i<name>] [-D<name>]" where "name" stands
 * for the string stored in the Descr field. Space for the help string is
 * allocated within.
 */
char *opt_helpstr(void);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER cmdopt.h ==== */

#endif	/* CMDOPT_H */
