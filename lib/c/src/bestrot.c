/* ==== FUNCTIONS bestrot.c ==== */

/* An implementation of the 3D point set alignment algorithm by
 * A. D. McLachlan. Reference:
 * McLachlan, A. D. (1979): J. Mol. Biol. 128: 49-79.
 * Replaces the buggy Kabsch rotation algorithm.
 */

/* ANSI C, 30. Apr. 1998. Andris Aszodi */

/* ---- HEADER ---- */

#include "bestrot.h"

/* ---- INCLUDE FILES ---- */

#include "siva.h"	    /* singular value decomposition */

/* ---- DEFINITIONS ---- */

#define DIM 3

/* ==== FUNCTIONS ==== */

/* center_vectors: calculates the weighted centroid of the 
 * set of 3-dimensional vectors X (Vno x 3) and
 * subtracts it from each of them,  thus centring the
 * set on the centroid. If Ctr==NULL, then a 3-long
 * array is allocated to store the centroid coordinates;
 * if Ctr!=NULL, then it is assumed to be large enough to
 * hold the coordinates.
 * If W==NULL, then uniform weighting is assumed.
 * Return value: Ctr, or NULL if Vno==0.
 */
double *center_vectors(double **X, double *Ctr, const double *W,
	 unsigned int Vno)
{
    register unsigned int i, j;
    register double Wsum;
    
    if (!Vno) return(NULL);
    
    /* allocate centroid vector if absent */
    if (Ctr==NULL)
	Ctr=(double *) calloc(DIM, sizeof(double));
    
    /* get the weighted centroid */
    if (W!=NULL)
    {
    	Wsum=0.0;
    	for(i=0; i<Vno; i++) Wsum+=W[i];
    }
    else Wsum=Vno;
    for (j=0; j<DIM; j++)
    {
	Ctr[j]=0.0;
	for (i=0; i<Vno; i++)
	    Ctr[j]+=(W!=NULL)? W[i]*X[i][j]: X[i][j];   /* suggested by WRT */
	Ctr[j]/=Wsum;
    }
    
    /* subtract Ctr from all vectors in X */
    for (i=0; i<Vno; i++)
	for (j=0; j<DIM; j++)
	    X[i][j]-=Ctr[j];
    
    return(Ctr);
}
/* END of center_vectors */

/* best_rot: finds the best rotation matrix that brings a set of
 * vectors X into another set Y. X, Y have Vno vectors (in rows), 
 * and both live in 3 dimensions (Vno x 3). W is a Vno-long
 * weight vector that can emphasise vector pairs. 
 * If W==NULL, then uniform weighting is assumed.
 * Transform is a 3x3 square matrix (allocated before call) that on
 * return contains the X->Y transformation. It is assumed that
 * X and Y were centered before the call.
 * NOTE: this routine cannot handle the degenerate cases when
 * the point sets are Dim<3-dimensional. (Might be implemented
 * later.) When this happens, a warning is printed to stderr
 * and -1.0 (a meaningless RMS value) is returned.
 * Return value: a weighted least-squares error function.
 */
double best_rot(double **X, double **Y, const double *W, 
		unsigned int Vno, Sqmat_ Transform)
{
    register unsigned int i, j, k, n, m;
    int Psign, Rank;
    double **H=NULL, *D=NULL, **K=NULL;
    Sqmat_ U=NULL;
    register double Err=0.0, Temp1, Temp2, Wsum=0.0;
    double Detu;
    
    /* set up the matrix to be SVD-d */
    U=alloc_sqmat(DIM);
    for (i=0; i<DIM; i++)
	for (j=0; j<DIM; j++)
	{
	    Temp1=0.0;
	    for (k=0; k<Vno; k++)
		Temp1+=(W!=NULL)? W[k]*X[k][i]*Y[k][j]: X[k][i]*Y[k][j];
	    U[i][j]=Temp1;
	}
    
    /* set up and perform SVD */
    siva_setup(DIM, DIM, &H, &D, &K);
    siva_decomp((const double **)U, DIM, DIM, H, D, K);
    
    /* check rank: do nothing if rank was lost */
    if ((Rank=rank_cond(D, DIM, SIVA_EPSILON, NULL))<3)
    {
	fprintf(stderr, "? best_rot(): Rank %d<%d\n", Rank, DIM);
	free_siva(DIM, DIM, H, D, K); free_matrix(U);
	return(-1.0);
    }
    
    /* get the determinant of U and store its sign in Psign */
    Psign=lu_decomp(U, DIM, NULL);
    Detu=lu_det(U, Psign, DIM);
    Psign=(Detu>0)? 1: -1;
    
    /* generate the transform matrix: here we explicitly
     * use DIM==3 because McLachlan does not say what to do
     * if Psign==-1 and DIM>3 and I don't know :-)
     */
    for (i=0; i<DIM; i++)
	for (j=0; j<DIM; j++)
	    Transform[i][j]=
		K[i][0]*H[j][0]+K[i][1]*H[j][1]+Psign*K[i][2]*H[j][2];
    
    /* evaluate the error function */
    for (n=0; n<Vno; n++)
    {
	Temp2=0.0;
	for (i=0; i<DIM; i++)
	{
	    Temp1=0.0;
	    for (j=0; j<DIM; j++)
		Temp1+=Transform[i][j]*X[n][j];
	    Temp1-=Y[n][i];
	    Temp2+=Temp1*Temp1;
	}
	if (W!=NULL)
	{
	    Err+=W[n]*Temp2;
	    Wsum+=W[n];
	}
	else Err+=Temp2;
    }
    Err/=(W!=NULL)? Wsum: Vno;
    
    /* cleanup */
    free_siva(DIM, DIM, H, D, K); free_matrix(U);
    return(sqrt(Err));
}
/* END of best_rot */

#undef DIM

/* ==== END OF FUNCTIONS bestrot.c ==== */
