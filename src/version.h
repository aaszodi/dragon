#ifndef VERSION_H
#define VERSION_H

/* ==== PROJECT DRAGON: HEADER version.h ==== */

/* Keeps track of the version string and the time/date
 * of last compilation/linking. This module is always recompiled
 * when a link is performed: the DRAGON_VERSION macro has to be passed
 * from the topmost Makefile. Compile as:-
 * $(CC) $(CFLAGS) -DDRAGON_VERSION='"$(VERSION)"' -c $(SRC)/version.c -o $(SRC)/version.o;
 * where $(VERSION) is the version macro string, $(SRC) is the C source directory.
 */

/* ANSI C, IRIX 5.3, 9. June 1995. Andras Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* version_string: returns a pointer to a (static) stringwhich contains
 * the program name, its version (coming from a macro passed through
 * the compiler command line) and the two ANSI C macros __DATE__ and __TIME__.
 */
char *version_string(void);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER version.h ==== */

#endif	/* VERSION_H */
