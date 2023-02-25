// ==== PROJECT DRAGON: METHODS Specgrad.c++ ====

/* Majorization algorithm employing the Spectral Gradient Method
 * (Glunt W, Hayden TL, Raydan M, J. Comput. Chem. 14:114-120 (1993)).
 */

// SGI C++ 7.1, IRIX 6.2, 9. May 1998. Andris Aszodi

// ---- HEADERS ----

#include <iostream.h>
#include <iomanip.h>
#ifdef __sgi
#include <ieeefp.h>
#endif
#include "Specgrad.h"

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this. Get around by the following macro
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#define fabsf fabs
#endif

// ---- Static initialisation ----

const double Specgrad_::SMALL=sqrt(DBL_MIN)/DBL_EPSILON;	// unsafe to do 1.0/SMALL

// ==== METHODS ====

/* weight(): sets up the calling object to work with a given
 * weight matrix W (with entries >=0.0).
 * Returns the size of the problem.
 */
int Specgrad_::weight(const Trimat_& W)
{
    N=W.rno();	    // save new problem size
    if (!N)
    {
	cerr<<"\n? Specgrad_::weight(): 0 matrix size\n";
	return(0);
    }
    Wgt=W; Wnorm=1.0;
    Distact.set_size(N); Bmat.set_size(N);
    Smat.set_size(N);
    make_smat();    // construct the S matrix from the weights
    
    return(N);
}
// END of init()

/* iterate(): performs the iteration on the point set Coords
 * (all vectors are assumed to have the same dimension). The coordinates
 * will be massaged towards the ideal UNsquared distances in Id. 
 * Eps is the relative precision (default 0.001).
 * Itno is the maximal number of iterations which contains
 * the actual iteration number on return. If no convergence was reached
 * then Itno is set to -Itno on return.
 * Return value: the "stress" (weighted dist difference).
 * Negative stress values indicate serious errors.
 */
float Specgrad_::iterate(const Trimat_& Id, Points_& Coords, 
	int& Itno, float Eps)
{
    // size checks
    if (Id.rno()<N)
    {
	cerr<<"\n? Specgrad_::iterate(): Ideal dist matrix dim too small ("
	    <<Id.rno()<<"<"<<N<<")\n";
	return(-1.0);
    }
    D=Coords.dim();
    if (!D)
    {
	cerr<<"\n? Specgrad_::iterate(): Dim mismatch within point set\n";
	return(-2.0);
    }
    if (Coords.active_len()<N)
    {
	cerr<<"\n? Specgrad_::iterate(): Too few points ("
	    <<Coords.active_len()<<"<"<<N<<")\n";
	return(-3.0);
    }
    Vector_ Ctr=Coords.centroid();
    Coords-=Ctr;
    norm_weights(Id);	// normalise weight matrix
    
    // set up the internal coordinate matrix and the negative gradient
    register unsigned int i, j;
    Xt.set_size(N, D); Xtbest.set_size(N, D);
    Negrad.set_size(N, D); Oldnegrad.set_size(N, D);
    for (i=0; i<N; i++)
	for (j=0; j<D; j++)
	    Xt[i][j]=Coords[i][j];
    
    // perform the iteration
    int Iter, Maxiter=Itno, Bkstep=0, Saveno=0;
    float Stress=stress(Id), Ostress=0.0, Bestress=Stress, Alpha=1.0;
    
    // bootstrap
    actual_dist();	// generate actual distances
    make_bmat(Id);	// get "B" matrix
    make_negrad();	// calc negative gradient
    Eps=fabsf(Eps);
    
    for (Iter=1; Stress>1e-6 && Iter<=Maxiter && Bkstep<=Maxiter; Iter++)
    {
	update_coords(Alpha);	// get new coordinates
	actual_dist();	// generate actual distances
	Ostress=Stress; Stress=stress(Id);
	if (Stress>=Ostress) { --Iter; ++Bkstep; } // didn't go up: do another step
	else
	{
	    Bkstep=0;	// reset backstep counter
	    if (Stress<Bestress)    // best so far, save
	    {
		Bestress=Stress;
		Xtbest=Xt; Saveno++;
	    }
	    if (fabsf(Stress-Ostress)<=Eps*Ostress)
		break;	// goes down and is good enough
	}
	
	Oldnegrad=Negrad;
	make_bmat(Id);
	make_negrad();	// here is the new gradient
	make_alpha(Alpha);  // and the new "stepsize"
    }

    // prepare results
    if (!Saveno)
    {
	// did not work
	cerr<<"\n? Specgrad_::iterate(Maxiter="<<Maxiter
	    <<", Eps="<<Eps<<"): No convergence\n";
	Itno=-Maxiter;	// negative iterno means failed convergence
    }
    else
    {
	// copy best coordinates back
	for (i=0; i<N; i++)
	    for (j=0; j<D; j++)
		Coords[i][j]=Xtbest[i][j];
	Itno=Iter;	// report back actual no. of iterations
    }
    Coords+=Ctr;	// shift to original centroid
    return(Bestress);
}
// END of iterate()

/* make_smat(): constructs the "S"-matrix from the
 * internal weight matrix Wgt. Protected
 */
void Specgrad_::make_smat()
{
    register unsigned int i, j;
    register double Temp;
    
    for (i=0; i<N; i++)
    {
	Smat[i][i]=0.0;
	for (j=0; j<i; j++)
	    Smat[i][j]=-Wgt[i][j];
    }
    for (i=0; i<N; i++)
    {
	Temp=0.0;
	for (j=0; j<N; j++) Temp+=(i>=j)? Smat[i][j]: Smat[j][i];
	Smat[i][i]=-Temp;
    }
}
// END of make_smat()

/* norm_weights(): norms the weights so that the weighted
 * squared sum of ideal distances in Id will equal 1.
 * The norm factor is stored in Wnorm. Protected
 */
void Specgrad_::norm_weights(const Trimat_& Id)
{
    Wgt*=Wnorm;	// reset previous state
    Smat*=Wnorm;
    
    register unsigned int i, j;
    register double Temp;
    Wnorm=0.0;
    for (i=0; i<N; i++)
	for (j=0; j<i; j++)
	{
	    if (Wgt[i][j]<=0.0) continue;
	    
	    Temp=Id[i][j];
	    Temp*=Temp;	// squared ideal dist
	    Temp*=Wgt[i][j];
	    Wnorm+=Temp;
	}
    if (Wnorm>SMALL)
    {
	Wgt/=Wnorm; Smat/=Wnorm;
    }
    else Wnorm=1.0;
}
// END of norm_weights()

/* actual_dist(): obtains the member matrix of UNsquared actual
 * distances (Distact) from the vectors in Xt. Protected
 */
void Specgrad_::actual_dist()
{
    // NOTE: the horrible spelt-out cycles are meant to save time
    register unsigned int i, j, k;
    register double Temp, Temp2;
    
    for (i=0; i<N; i++)
    {
	Distact[i][i]=0.0;
	for (j=0; j<i; j++)
	{
	    Temp=0.0;
	    for (k=0; k<D; k++)
	    {
		Temp2=Xt[i][k]-Xt[j][k];
		Temp+=Temp2*Temp2;
	    }
	    Distact[i][j]=sqrt(Temp);
	}
    }
}
// END of actual_dist()

/* make_bmat(): constructs the "B" matrix from the weights and the
 * ideal and actual (unsquared) distances. Protected
 */
void Specgrad_::make_bmat(const Trimat_& Distid)
{
    register unsigned int i, j;
    register double Temp, Temp2;
    
    for (i=0; i<N; i++)
    {
	Temp=0.0;
	// rows
	for (j=0; j<i; j++)
	{
	    Temp2=(Wgt[i][j]<=SMALL || Distact[i][j]<=SMALL)? 
		0.0: -Wgt[i][j]*Distid[i][j]/Distact[i][j];
	    Bmat[i][j]=Temp2;
	    Temp+=Temp2;
	}
	// cols
	for (j=i+1; j<N; j++)
	{
	    Temp2=(Wgt[j][i]<=SMALL || Distact[j][i]<=SMALL)? 
		0.0: -Wgt[j][i]*Distid[j][i]/Distact[j][i];
	    Bmat[j][i]=Temp2;
	    Temp+=Temp2;
	}
	Bmat[i][i]=-Temp;   // diagonal sum
    }
}
// END of make_bmat()

/* stress(): calculates the "stress" value, i.e. the weighted
 * difference between the ideal and actual distances. Protected
 */
float Specgrad_::stress(const Trimat_& Distid) const
{
    register unsigned int i, j;
    register float Stress=0.0, Temp;
    
    for (i=0; i<N; i++)
	for (j=0; j<i; j++)
	{
	    if (Wgt[i][j]<=0.0) continue;
	    Temp=float(Distid[i][j])-float(Distact[i][j]);
	    Stress+=float(Wgt[i][j])*Temp*Temp;
	}
    return(Stress);
}
// END of stress()

/* make_negrad(): makes the negative gradient of the stress function. Protected */
void Specgrad_::make_negrad()
{
    // WARNING: Bmat is modified!
    Bmat-=Smat;

    register unsigned int i, j, k;
    register double Temp;
    
    for (i=0; i<N; i++)
	for (j=0; j<D; j++)
	{
	    Temp=0.0;
	    for (k=0; k<N; k++)
		Temp+=Bmat(i, k)*Xt[k][j];
	    Negrad[i][j]=2.0*Temp;
	}
}
// END of make_negrad()

/* update_coords(): updates the coordinates with the negative gradient,
 * using Alpha as the "stepsize". Protected
 */
void Specgrad_::update_coords(float Alpha)
{
    if (!finite(Alpha) || fabs(Alpha)<SMALL)
    {
	cerr<<"\n? Specgrad_::update_coords("<<Alpha<<")"<<endl;
	return;
    }
    Alpha=1.0/Alpha;
    
    register unsigned int i, j;
    for (i=0; i<N; i++)
	for (j=0; j<D; j++)
	    Xt[i][j]+=Alpha*Negrad[i][j];
}
// END of update_coords()

/* make_alpha(): calculates the new stepsize Alpha. Protected */
void Specgrad_::make_alpha(float& Alpha) const
{
    register unsigned int i, j;
    register double Num=0.0, Denom=0.0, Temp;
    
    // generate inner products
    for (i=0; i<N; i++)
	for (j=0; j<D; j++)
	{
	    Temp=Oldnegrad[i][j];
	    Num+=Negrad[i][j]*Temp;
	    Denom+=Temp*Temp;
	}
    if (Denom>SMALL) Alpha*=(1.0-Num/Denom);
}
// END of make_alpha()

// ==== END OF METHODS Specgrad.c++ ====
