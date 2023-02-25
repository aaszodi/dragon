// ==== PROJECT DRAGON: FUNCTIONS Hmom.c++ ====

/* Calculates the hydrophobic moments of clusters and
 * rotates them in Euclidean space so that the moments point
 * towards the common centroid.
 */

// SGI C++ 7.1, IRIX 6.2, 9. May 1998. Andris Aszodi 

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

// ---- UTILITY HEADERS ----

#include "Vector.h"
#include "Sqmat.h"
#include "Svd.h"

// ---- MODULE HEADERS ----

#include "Hmom.h"
#include "Fakebeta.h"

// ---- DEFINITIONS ----

#ifdef M_PI
#define TWO_PI (2.0*M_PI)
#else
#define TWO_PI 6.2831853
#endif

// ---- PROTOTYPES ----

static int rot_ndim(const Vector_& P, const Vector_& Q, Sqmat_& R);

// ==== FUNCTIONS ====

/* hmom_clurot(): rotates all clusters in Pieces so that their
 * hydrophobic moments point towards the overall centroid.
 * The phobicity info is in Polymer, the Euclidean coordinates
 * are in Xyz. No action is taken for single-cluster sets.
 */
void hmom_clurot(const Pieces_& Pieces, const Polymer_& Polymer, 
	Points_& Xyz)
{
    static Points_ Beta;   // fake beta positions
    static Vector_ Ctr, Ab, Hmom;
    static Sqmat_ Rot;
    
    unsigned int Ptno=Xyz.len(), Cluno=Pieces.clu_no();
    if (Cluno<=1) return;
    
    Bits_ Oldmask=Xyz.mask(true);   // switch all points on
    unsigned int Dim=Xyz.dim();
    
    if (!Dim)
    {
	cerr<<"\n? hmom_clurot(): Dim mismatch within coords\n";
	return;
    }
    
    // get fake beta positions
    Beta.len_dim(Ptno, Dim);
    Fakebeta_::beta_xyz(Xyz, Polymer, Beta);
    
    // store moment vectors in Beta
    Ab.dim(Dim);
    register unsigned int i;
    
    for (i=0; i<Ptno; i++)
    {
	if (i==0 || i==Ptno-1)	// no fake beta on these (N/C termini
	{
	    Beta[i].set_values();   // just zero it
	    continue;
	}
	
	Ab=Beta[i]; Ab-=Xyz[i];
	Ab.vec_norm();
	Ab*=Polymer.phob(i-1);	// shift because of N/C termini
	Beta[i]=Ab;
    }
    
    // adjust all clusters one by one
    register unsigned int ci;
    
    Ctr.dim(Dim); 
    Rot.set_size(Dim); Hmom.dim(Dim); 
    for (ci=0; ci<Cluno; ci++)
    {
	if (Pieces.clus(ci).on_no()<=1)		// 1-size cluster: do nothing
	    continue;
	
	// check cluster type and skip coils: 4-Apr-1997.
	Pieces_::Clutype_ Clutyp=Pieces.clu_type(ci);
	if (Clutyp==Pieces_::COIL)
	    continue;
	
	Xyz.mask(Pieces.clus(ci));
	Beta.mask(Pieces.clus(ci));

	// sum the moment vectors in Beta for current cluster
	Hmom.set_values();  // zero
	for (i=0; i<Beta.active_len(); i++)
	    Hmom+=Beta[i];

	Ctr=Xyz.centroid(); Ctr*=-1.0;	// point towards overall centroid
	if (!rot_ndim(Hmom, Ctr, Rot))	// get DAMPED rotation matrix
	    continue;    // skip rest on error
	Xyz+=Ctr;   // center cluster: observe sign change!
	Xyz*=Rot;   // rotate so that Hmom points to overall centroid 
	Xyz-=Ctr;   // move back to original place
    }
    Xyz.mask(Oldmask);
}
// END of hmom_clurot()

/* rot_ndim(): constructs a rotation matrix which
 * rotates the vector P into another vector Q, leaving the
 * subspace orthogonal to the 2-dimensional subspace spanned
 * by P and Q intact. Returns the square matrix R (size adjusted
 * within if necessary) that rotates P onto Q. Rotation means that
 * the direction of P is changed into the direction of Q, 
 * so |P|==|Q| is not necessary. In case of errors
 * (dim mismatch, collinearity etc.) R is a unit matrix on return.
 * NOTE: from Version 4.4.2 on, this routine contains a damping
 * and therefore the rotations are not exact!!! See comments inside
 * Return value: 1 if OK, 0 on error.
 */
static int rot_ndim(const Vector_& P, const Vector_& Q, Sqmat_& R)
{
    static Sqmat_ B, U;
    static Svd_ Svd;
    
    unsigned int N=P.dim();
    B.set_size(N); B.diag_matrix();   // NxN unit matrix
    
    // paranoid section
    if (N<2)
    {
	cerr<<"\n? rot_ndim(): "<<N<<"-dim case\n";
	return(0);
    }
    
    if (N!=Q.dim())
    {
	cerr<<"\n? rot_ndim(): P, Q dim mismatch\n";
	return(0);
    }
    
    // resize output and init to diag matrix
    R.set_size(N);
    R.diag_matrix();
    
    // make unit vectors, check for null-vectors
    Vector_ Pu(P);
    double Plen=Pu.vec_norm();
    if (Plen==0.0)
    {
	cerr<<"\n? rot_ndim(): P is nullvector\n";
	return(0);
    }
    Vector_ Qu(Q);
    double Qlen=Qu.vec_norm();
    if (Qlen==0.0)
    {
	cerr<<"\n? rot_ndim(): Q is nullvector\n";
	return(0);
    }
    
    double Cosfi=Pu*Qu;
    if (1.0-SVD_EPSILON<Cosfi && Cosfi<1.0+SVD_EPSILON)	// Pu==Qu
    {
	cerr<<"\n? rot_ndim(): P, Q are collinear, <P|Q>="<<Cosfi<<endl;
	return(1);
    }
    
    if (-1.0-SVD_EPSILON<Cosfi && Cosfi<-1.0+SVD_EPSILON)  // Pu==-Qu
    {
	cerr<<"\n? rot_ndim(): (PQ) angle is 180 degrees, <P|Q>="<<Cosfi<<endl;
	return(1);
    }
    
    /* Get the two orthonormal base vectors which span the
     * P:Q plane. This is, in fact, a little Gram/Schmidt
     * orthogonalisation
     */
    Vector_ B1(Pu);  // already unit vector
    Vector_ B2=Qu-B1*Cosfi;
    B2.vec_norm();
    
    /* Fill up the B matrix with the "orthocomplements" of the
     * original base vectors, i.e. decompose each column in B
     * (which is the unit matrix here) into a "parallel part"
     * that lies in the P:Q plane and a complementary "orthogonal part"
     * and replace each column with the orthogonal part. The
     * resulting matrix should have a rank of N-2.
     */
    Vector_ Ev(N);
    register unsigned int i, k;
    for (i=0; i<N; i++)
    {
	Ev=B.col(i);	// get i-th base vector
	Ev-=(B1*(Ev*B1)+B2*(Ev*B2));	// subtract P:Q plane component
	if (Ev.vec_norm()<SVD_EPSILON)
	{
	    cerr<<"\n? rot_ndim(): B column ["<<i<<"] is almost 0-vector\n";
	    // do something here
	}
	B.col(Ev, i);	// put back orthogonal part
    }
    // do an SVD which constructs an orthonormal basis from this
    Svd.set_size(N, N);
    if (Svd.make_decomp(B))
    {
	cerr<<"\n? rot_ndim(): SVD decomposition error\n";
	return(0);
    }

    /* scan for the two smallest W values: their indices are the
     * col indices to be replaced by the base vectors B1, B2.
     */
    register double Wi, Sm=HUGE_VAL, Sm2=HUGE_VAL;
    register unsigned int Si, Si2;
    for (i=0; i<N; i++)
    {
	Wi=Svd.w()[i];
	if (Wi<Sm)	// new smallest
	{
	    Sm2=Sm; Sm=Wi;
	    Si2=Si; Si=i;
	}
	else if (Wi<Sm2)    // new second smallest
	{
	    Sm2=Wi; Si2=i;
	}
    }

    /* copy the N-2 largest-W cols from u() into U and
     * pop the B1,B2 base vectors into the LAST two columns of U
     * (the first N-2 cols will span P:Q's orthocomplement)
     */
    U.set_size(N);
    for (i=k=0; i<N; i++)
    {
	if (i==Si || i==Si2) continue;	// don't copy the two smallest
	Ev=Svd.u().col(i);
	U.col(Ev, k++);
    }
    U.col(B1, N-2); U.col(B2, N-1);	// last two
    /* Get the 2D rotation matrix. The rotation angle Phi will be
     * damped so that as it approaches Pi, the damped angle approaches
     * Pi/2. IMPORTANT NOTE: if you ever plan to use this routine
     * outside DRAGON-IV for exact rotations,  excise the following
     * part of code!
     */
    double Phi=acos(Cosfi);
    Phi=Phi*(TWO_PI-Phi)/TWO_PI;
    Cosfi=cos(Phi);
    /* <<<< END OF ROTATION ANGLE DAMPING >>>> */
    
    double Sinfi=sqrt(1.0-Cosfi*Cosfi);

    /* get the overall matrix: R=U*R2*U' where R2 is a unit matrix
     * except for the lower right corner which is an "angle Phi
     * 2D rotation block",  i. e.
     * [ cos(Phi) -sin(Phi)]
     * [ sin(Phi) cos(Phi) ]
     * The multiplication is expanded for speed (?) and obfuscation.
     */
    register unsigned int m, n, p;
    register double Temp;
    register double *Um, *Un;	// U row ptrs (m-th and n-th) for speed
    for (m=0; m<N; m++)
    {
	Um=U[m];
	for (n=0; n<N; n++)
	{
	    Un=U[n];
	    Temp=0.0;
	    for (p=0; p<=N-3; p++) Temp+=(Um[p]*Un[p]);
	    Temp+=Um[N-2]*(Cosfi*Un[N-2]-Sinfi*Un[N-1])+
		Um[N-1]*(Sinfi*Un[N-2]+Cosfi*Un[N-1]);
	    R[m][n]=Temp;
	}
    }
    return(1);
}
// END of rot_ndim()

// ==== END OF FUNCTIONS Hmom.c++ ====
