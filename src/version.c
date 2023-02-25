/* ==== PROJECT DRAGON: FUNCTIONS version.c ==== */

/* Keeps track of the version string and the time/date
 * of last compilation/linking. This module is always recompiled
 * when a link is performed: the DRAGON_VERSION macro has to be passed
 * from the topmost Makefile. Compile as:-
 * $(CC) $(CFLAGS) -DDRAGON_VERSION='"$(VERSION)"' -c $(SRC)/version.c -o $(SRC)/version.o;
 * where $(VERSION) is the version macro string, $(SRC) is the C source directory.
 */

/* ANSI C, IRIX 5.3, 9. June 1995. Andras Aszodi */

#ifndef DRAGON_VERSION
#error Please specify the DRAGON_VERSION macro
#endif

/* ---- MODULE HEADER ---- */

#include "version.h"

/* ==== FUNCTIONS ==== */

/* version_string: returns a pointer to a (static) string which contains
 * the program name, its version (coming from a macro passed through
 * the compiler command line) and the two ANSI C macros __DATE__ and __TIME__.
 */
char *version_string(void)
{
    static char Vstr[80];
    static int Init=1;
    
    if (Init)
    {
	sprintf(Vstr, "DRAGON %s [%s, %s]", DRAGON_VERSION, __DATE__, __TIME__);
	Init=0;
    }
    return(Vstr);
}
/* END of version_string */

/* ==== END OF FUNCTIONS version.c ==== */
