// ==== PROJECT DRAGON: METHODS Sstrbase.c++ ====

/* Base class for handling secondary structures: H-bond topology
 * and ideal geometry. Chain topology comes from the Segment module.
 */

// SGI C++ 4.0, IRIX 5.3, 2. Oct. 1995. (C) Andras Aszodi

// ---- STANDARD HEADERS ----

#include <math.h>
#include <string.h>

// ---- UTILITY HEADERS ----

#include "Vector.h"
#include "Sqmat.h"
#include "Hirot.h"

// ---- MODULE HEADER ----

#include "Sstrbase.h"

// ---- DEFINITIONS ----

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

// ==== Sstrbase_ METHODS ====

/* make_helix(): constructs an ideal helix in a point set Hel which
 * must be properly masked and the active region should be 3-dimensional.
 * The helix parameters are supplied by Radius, Pitch and Turn and
 * the structure will be grown so that the N->C direction corresponds
 * to the positive direction on the X axis. Phasing (default 1)
 * determines whether the first point is on +Y (>0) or -Y (when <=0).
 * Return value: the length of the helix built or 0 on error. Protected static
 */
unsigned int Sstrbase_::make_helix(Points_& Hel, double Radius, 
	double Pitch, double Turn, int Phasing)
{
    if (Hel.dim()!=3)
    {
	cerr<<"\n? Sstrbase_::make_helix(): Sorry,  not 3D\n";
	return(0);
    }
    
    unsigned int L=Hel.active_len();
    if (!L) return(0);	// totally inactive
    
    // flip structure around X-axis (makes sense for beta-strands only)
    if (Phasing<=0) Radius*=(-1);
    
    for (register unsigned int i=0; i<L; i++)
    {
	Hel[i][0]=i*Pitch;  // X, Y, Z coords
	Hel[i][1]=Radius*cos(i*Turn);
	Hel[i][2]=Radius*sin(i*Turn);
    }
    return(L);
}
// END of make_helix()

/* pos4_angle: given the Cartesian coordinates of 4 points (P[1-4]),
 * the torsion angle defined by them (along 2-3) is returned. The
 * value should be between -Pi..+Pi with the sign indicating the
 * usual handedness convention. -2*Pi is returned if any 3 of the 
 * points are colinear. Protected static
 */
double Sstrbase_::pos4_angle(const Vector_& P1, const Vector_& P2,
	const Vector_& P3, const Vector_& P4)
{
    // get normal vectors of planes (123) and (234)
    Vector_ V2=P3-P2;   // torsion is taken along this
    Vector_ W1=cross_prod(P2-P1,V2);
    Vector_ W2=cross_prod(V2,P4-P3);

    // colinearity check: return -2Pi if fails
    register double W1len, W2len;

    if (0.0==(W1len=W1.vec_len()) ||
	0.0==(W2len=W2.vec_len()))
	    return(-2.0*M_PI);

    // get torsion angle (i.e. angle between planes (123),(234))
    register double Costheta;
    
    Costheta=(W1*W2)/(W1len*W2len);
    Costheta=acos(Costheta);

    // get sign of theta (handedness)
    return(((V2*cross_prod(W1, W2))>=0.0)? Costheta: -Costheta);
}
// END of pos4_angle()

/* <<: calls a virtual output function. Details in derived classes */
ostream& operator<<(ostream& Out, const Sstrbase_& S)
{
    S.write_to(Out);
    return(Out);
}
// END of <<

// ==== END OF METHODS Sstrbase.c++ ====
