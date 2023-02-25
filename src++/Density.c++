// ==== PROJECT DRAGON: FUNCTIONS Density.c++ ====

/* Adjusts molecular density. */

// SGI C++ 4.0, IRIX 5.3, 4. July 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

// ---- MODULE HEADER ----

#include "Density.h"

// ---- UTILITY HEADERS ----

#include "Sqmat.h"
#include "Ql.h"

// ---- PROTOTYPES ----

static void update_coords(double Densfact, const Pieces_& Pieces, 
	Points_& Xyz);
static int dbl_cmp(const void *X, const void *Y);

// ==== FUNCTIONS ====

// ---- Distance-based adjustments ----

/* scale_distdens(): uses Willie's "moment scaling" to adjust the
 * density of the distance matrix to the density expected for a
 * 3D spherical protein. He has calculated the first and second moments
 * of the distance distribution between any pairs of points within
 * a solid sphere of radius Rmax. Here the moments are found and a
 * scale factor is calculated as the average of the ratios of the
 * expected moments to the observed ones.
 * Return value: the scaling factor.
 */
double scale_distdens(Trimat_& Dist, double Rmax)
{
    register double Dij, Davg, Davg2, Dexp, Dexp2;
    double Densfactor;
    register unsigned int i, j, d, No, Ptno=Dist.rno();

    // init
    No=(Ptno*(Ptno-1))/2;     // excludes the main diagonal
    if (!No) return(1.0);   // chain too short

    Dexp=36.0*Rmax/35.0; Dexp2=1.2*Rmax*Rmax; // avg and SD
    Davg=Davg2=0.0;

    /* average distances and squared distances: because Dist[][]
     * contains the squared values, the square root is taken for
     * the distance average. Main diag skipped
     */
    for (d=1; d<Ptno; d++)
        for (i=d; i<Ptno; i++)
        {
            j=i-d;
            Dij=Dist[i][j];
            Davg2+=Dij;
            Davg+=sqrt(Dij);
        }
    Davg/=No; Davg2/=No;

    /* the factor is the average of the ideal/observed
     * ratios of the 1st and 2nd moments. 
     */
    Densfactor=((Dexp/Davg)+sqrt(Dexp2/Davg2))/2.0;

    /* adjust the density of Dist: square Densfactor first */
    Densfactor=Densfactor*Densfactor;
    for (d=2; d<Ptno; d++)               // all except main and first off-diag
        for (i=d; i<Ptno; i++)
        {
            j=i-d;
            Dist[i][j]*=Densfactor;
        }
    return(Densfactor);
}
// END of scale_distdens()

/* proj_dens(): projection generally shrinks the coordinates. Here
 * a new distance matrix is computed from the Xyz[][] coords, and
 * compared to the old Dist[][]. The adjustment factor Fact is
 * chosen so that SUM(Dist[i][j]-Fact*Newdist[i][j])^2 be minimal
 * (very simple linear regression). 
 * If the point set is composed of multiple clusters, then the
 * clusters are moved away from each other as rigid bodies:
 * otherwise the whole assembly is "blown up". The distance
 * matrix is not changed.
 * Return value: sqrt(Fact).
 */
double proj_dens(const Trimat_& Dist, const Pieces_& Pieces, Points_& Xyz)
{
    static Trimat_ Newdist; // will be resized every time when necessary

    double Fact, Sdx, Sx2, Xij;
    register unsigned int i, j,Ptno=Dist.rno();

    // get new dist matrix from Xyz
    Xyz.dist_mat2(Newdist);

    // linear regression: Dist[i][j]=Fact*Newdist[i][j] if possible
    Sdx=Sx2=0.0;
    for (i=1; i<Ptno; i++)
        for (j=0; j<i; j++)
        {
            Xij=Newdist[i][j];
            Sdx+=sqrt(Dist[i][j]*Xij);
            Sx2+=Xij;
        }
    Fact=Sdx/Sx2;

    // update coordinates
    update_coords(Fact, Pieces, Xyz);
    return(Fact);
}
// END of proj_dens()

// ---- 3D Euclidean adjustments ----

/* ellips_dens(): adjusts the density of the 3D Euclidean point set Xyz
 * to match the expected density Expdens. Fits an ellipsoid to the points
 * so that it contains 90% of them,  works out an adjustment factor and
 * updates the points. Single-cluster point sets are isotropically adjusted, 
 * multicluster sets are updated by moving the clusters as rigid bodies.
 * Return value: the adjustment factor. If the dimension is not 3D, 
 * then no action is taken and 0.0 is returned.
 */
double ellips_dens(double Expdens, const Pieces_& Pieces, Points_& Xyz)
{
    if (Xyz.dim()!=3) return(0.0);  // no action
    
    // these used to be static but we can afford the allocation
    // to get less memory leak... 27-Apr-1998.
    Trimat_ Moment(3);
    Vector_ Evals(3);
    Sqmat_ Evec(3);
    double *Ellips=NULL;
    unsigned int Ptno=Xyz.len();
    
    Ellips=new double [Ptno];
    
    register unsigned int i, j, k;
    register double Temp;
    
    // fill up moment matrix
    Xyz.mask(true);	// switch on
    for (i=0; i<3; i++)
	for (j=0; j<=i; j++)
	{
	    Temp=0.0;
	    for (k=0; k<Ptno; k++) Temp+=Xyz[k][i]*Xyz[k][j];
	    Moment[i][j]=Temp;
	}
    
    // diagonalise
    if (eigen_ql(Moment, Evals, Evec))
    { free(Ellips); return(0.0);}	// iter limit exceeded (something nasty)
    
    // rotate to new coord system and substitute into 3D ellipsoid equation
    double Coord123;
    for (i=0; i<Ptno; i++)
    {
	Ellips[i]=0.0;
	for (j=0; j<3; j++)
	{
	    Coord123=Xyz[i]*Evec.row(j);
	    Ellips[i]+=(Coord123*Coord123/Evals[j]);
	}
    }
    
    // sort to get 90% containment
    qsort(Ellips, Ptno, sizeof(double), dbl_cmp);
    double Ellfact=Ellips[int(0.9*Ptno)-1];
    
    // get ellipsoid volume and density
    const double PI_43=4.1887902;	// 4*Pi/3
    double Vol, Density, Densfact;
    
    Vol=PI_43*sqrt(Ellfact*Ellfact*Ellfact*Evals[0]*Evals[1]*Evals[2]);
    Density=Ptno/Vol;	// residues per volume
    Densfact=Expdens/Density;	// cubic density factor
    Densfact=pow(Densfact, -1.0/3.0);	// linear factor
    
    // adjust density
    update_coords(Densfact, Pieces, Xyz);
    free(Ellips);
    return(Densfact);
}
// END of ellips_dens()

// ---- Auxiliaries ----

/* update_coords(): updates the Euclidean coordinates in Xyz with the
 * density adjustment factor Densfact. The cluster layout is in Pieces.
 * For single-cluster point sets, the whole set is inflated/deflated
 * as in 3.x; for multiclusters, the clusters are translated to/from the
 * common centroid as rigid bodies. Full activation is reset for Xyz
 * on return.
 */
static void update_coords(double Densfact, const Pieces_& Pieces, 
	Points_& Xyz)
{
    if (Pieces.clu_no()<=1)
	Xyz*=Densfact;	// do on all points (isotropic adjustment)
    else
    {
	Vector_ Ctr;
	for (register unsigned int i=0; i<Pieces.clu_no(); i++)
	{
	    Xyz.mask(Pieces.clus(i));	// do it cluster by cluster
	    Ctr=Xyz.centroid();		// get cluster centroid
	    Xyz+=((Densfact-1.0)*Ctr);	// move the cluster by Densfact
	}
	Xyz.mask(true);			// reset full activation
    }
}

/* dbl_cmp(): auxiliary function for the qsort() in ellips_dens(). */
static int dbl_cmp(const void *X, const void *Y)
{
    register double Temp;

    Temp= *((double *)X)- *((double *)Y);
    return((Temp>0.0)? 1: ((Temp<0.0)? -1: 0));
}
/* END of dbl_cmp() */

// ==== END OF FUNCTIONS Density.c++ ====


