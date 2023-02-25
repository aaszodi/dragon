#ifndef PORTRANDOM_HEADER
#define PORTRANDOM_HEADER

/* ==== HEADER portrandom.h ==== */

/* Portable random number generator. Based on "ran1" in
 * Numerical Recipes. Slightly altered to avoid the clumsy
 * initialisation.
 * Reference:
 * Numerical Recipes,  Second Edition, 1992 (ver. 2.02)
 * Chapter 7,  Page 280.
 */

/* ANSI C, IRIX 5.2, 27. Feb. 1996. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* init_portrand: initialises the portable random number generator.
 * If Seed==0, then Seed==1 is assumed.
 */
void init_portrand(long Seed);

/* port_rand: returns a non-negative long pseudo-random number.
 * Maximum number of sequential calls is around 10^8.
 */
long port_rand(void);

/* port_random: the portable random number generator itself.
 * Returns a pseudo-random number in the interval (0.0 .. 1.0).
 * Maximum number of sequential calls is around 10^8.
 */
double port_random(void);

/* portrandom_gauss: returns normally distributed random numbers with
 * zero mean and unit variance. Based on the Box/Muller method
 * as described in Numerical Recipes. 
 */
double portrandom_gauss(void);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER portrandom.h ==== */

#endif	/* PORTRANDOM_HEADER */
