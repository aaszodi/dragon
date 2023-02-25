/* ==== FUNCTIONS portrandom.c ==== */

/* Portable random number generator. Based on "ran1" in
 * Numerical Recipes. Slightly altered to avoid the clumsy
 * initialisation.
 * Reference:
 * Numerical Recipes,  Second Edition, 1992 (ver. 2.02)
 * Chapter 7,  Page 280.
 */

/* ANSI C, IRIX 4.0.5, 27. Feb. 1996. Andris Aszodi */

/* ---- HEADER ---- */

#include "portrandom.h"

/* ---- PRIVATE CONSTANTS AND VARIABLES ---- */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS (2.2e-15)
#define RNMX (1.0-EPS)

static long iy=0L;
static long iv[NTAB];
static long idum=-1L;

static char Spare=0;	/* Gauss variables */
static double Spval;

/* ==== FUNCTIONS ==== */

/* init_portrand: initialises the portable random number generator.
 * If Seed==0, then Seed==1 is assumed.
 */
void init_portrand(long Seed)
{
    register int j;
    register long k;
    
    idum=(!Seed)? 1L: Seed; /* do not init with 0 */
    if (idum<0) idum*=-1L;
    
    /* fill up the shuffle table */
    for (j=NTAB+7; j>=0; j--)
    {
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum<0) idum+=IM;
	if (j<NTAB) iv[j]=idum;
    }
    iy=iv[0];
    
    Spare=0;	/* unset previous value from portrandom_gauss() */
}
/* END of init_portrand */

/* port_rand: returns a non-negative long pseudo-random number.
 * Maximum number of sequential calls is around 10^8.
 */
long port_rand(void)
{
    register int j;
    register long k;
    register double Temp;
    
    if (idum<0 || !iy) init_portrand(1);    /* init if necessary */
    
    /* pick a random position (j) from the shuffle table iv[] */
    k=idum/IQ;
    idum=IA*(idum-k*IQ)-IR*k;
    if (idum<0) idum+=IM;
    j=iy/NDIV;
    iy=iv[j]; iv[j]=idum;   /* refill the shuffle table */
    return(iy);
}
/* END of port_rand */

/* port_random: the portable random number generator itself.
 * Returns a pseudo-random number in the interval (0.0 .. 1.0).
 * Maximum number of sequential calls is around 10^8.
 */
double port_random(void)
{
    register int j;
    register long k;
    register double Temp;
    
    /* this is the same as port_rand() but "inlined" */
    if (idum<0 || !iy) init_portrand(1);    /* init if necessary */
    
    /* pick a random position (j) from the shuffle table iv[] */
    k=idum/IQ;
    idum=IA*(idum-k*IQ)-IR*k;
    if (idum<0) idum+=IM;
    j=iy/NDIV;
    iy=iv[j]; iv[j]=idum;   /* refill the shuffle table */
    
    /* scale to 0.0 .. 1.0: return value should not be the endpoint */
    return(((Temp=AM*iy)>RNMX)? RNMX: Temp);
}
/* END of port_random */

/* portrandom_gauss: returns normally distributed random numbers with
 * zero mean and unit variance. Based on the Box/Muller method
 * as described in Numerical Recipes. 
 */
/* Machine precision constant comes from DBL_EPSILON in <float.h>
 * in pure ANSI C header environments (not on Sun4/SunOS4.x)
 */
#ifndef DBL_EPSILON
#define EPSILON 2.2e-15
#else
#define EPSILON DBL_EPSILON
#endif
double portrandom_gauss(void)
{
    register double Fac, R, V1, V2;
    
    if (Spare)	/* we had a value from a previous call: return it */
    {
	Spare=0;
	return(Spval);
    }
    else    /* make two values, save one, return the other */
    {
	do  /* get two random numbers within the unit circle */
	{
	    V1=2.0*port_random()-1.0; V2=2.0*port_random()-1.0;
	    R=V1*V1+V2*V2;
	}
	while (R>=1.0 || R<=EPSILON);
	Fac=sqrt(-2.0*log(R)/R);
	Spval=V1*Fac;	/* spare value to be returned next time */
	Spare=1;
	return(V2*Fac);
    }
}
#undef EPSILON
/* END of random_gauss */

#undef IA 
#undef IM 
#undef AM 
#undef IQ 
#undef IR 
#undef NTAB
#undef NDIV 
#undef EPS
#undef RNMX

/* ==== END OF FUNCTIONS portrandom.c ==== */
