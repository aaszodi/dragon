/* ==== PROJECT DRAGON: FUNCTION matherr.c ==== */

/* System V math exception handling. On SGI-s, link with the
 * -lmx option (extended math library). There is no header file:
 * the prototype for matherr() is in <math.h>.
 */

/* ANSI C+extensions, IRIX 5.3, 11. July 1996. Andris */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <float.h>
#define __EXTENSIONS__
#include <math.h>

/* matherr(): under System V, this function, when defined,
 * is invoked at certain floating-point exceptions. There is
 * no way at present to intercept divisions-by-zero, but the
 * result shows up sooner or later as a NaN.
 * DOMAIN errors are trapped and the program is aborted unless
 * the error occurred because of a negative argument to sqrt()
 * or sqrtf() in which case the absolute value of the argument
 * is used. UNDERFLOW and OVERFLOW are also trapped and the return
 * values are set to 0.0 or DBL_MAX, respectively. All other errors
 * are treated by the default mechanism in the library.
 */
int matherr(register struct exception *Ex)
{
    switch(Ex->type)
    {
	case DOMAIN:
	if (isnan(Ex->arg1) || isnan(Ex->arg2))	/* abort on NaN-s */
	{
	    fprintf(stderr, "\n! DOMAIN NaN argument(s) for %s, dying\n", Ex->name);
	    kill(SIGFPE, getpid()); /* kill myself */
	}
	else if (!strncmp(Ex->name, "sqrt", 4))  /* replace arg with abs val */
	{
	    fprintf(stderr, "\n? DOMAIN fp exception: %s(%f), abs val used\n", Ex->name, Ex->arg1);
	    Ex->retval=sqrt(-(Ex->arg1));
	}
	else
	{
	    fprintf(stderr, "\n! DOMAIN fp exception: op=%s, arg1=%f, arg2=%f, dying\n", 
		Ex->name, Ex->arg1, Ex->arg2);
	    kill(SIGFPE, getpid());
	}
	return(1);  /* no default handling */
	
	case UNDERFLOW:
	fprintf(stderr, "\n? UNDERFLOW fp exception: op=%s, arg1=%f, arg2=%f, result 0.0\n", 
	    Ex->name, Ex->arg1, Ex->arg2);
	Ex->retval=0.0;
	return(1);

	case OVERFLOW:
	fprintf(stderr, "\n? OVERFLOW fp exception: op=%s, arg1=%f, arg2=%f, result DBL_MAX\n", 
	    Ex->name, Ex->arg1, Ex->arg2);
	Ex->retval=DBL_MAX;
	return(1);
	
	default: return(0); /* handle as usual */
    }
}
/* END of matherr() */

/* ==== END OF FUNCTION matherr.c ==== */

