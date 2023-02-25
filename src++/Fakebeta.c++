// ==== PROJECT DRAGON: METHODS Fakebeta.c++ ====

/* Distances between the C-alphas on the backbone and the
 * fake C-beta atoms (representing the side-chains).
 * The C-beta positions are determined by the backbone.
 */

// SGI C++ 4.0, IRIX 5.3, 28. Nov. 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <iostream.h>
#include <iomanip.h>

// ---- CLASS HEADER ----

#include "Fakebeta.h"

// ---- FLOAT/DOUBLE ISSUES ----

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this. Get around by the following macro
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#define fabsf fabs
#endif

// ==== Fakebeta_ METHODS ====

// ---- Inline auxiliaries ----

/* get_dist: all distances in this scheme can be calculated 
 * from three other distances (D1..D3) and L (a Lambda).
 * That's all what this tiny function does. If L==1.0 then
 * the corresponding C-alpha and C-beta are at the same place;
 * in this case, D1 is returned. Static private
 */
inline
float Fakebeta_::get_dist(float D1, float D2, float D3, float L)
{
    float L1=1.0-L;
    return((L1==0.0)? D1: (D1-L*D2+L*L1*D3)/L1);
}
// END of get_dist()

// ---- Distance update ----

/* update(): updates the CA:CB and CB:CB distance matrices using the
 * CA:CA matrix Dista and the prescribed CA(i):CB(i) distances from
 * Polymer_ . The matrices within may be resized if necessary.
 * Return value: the new size.
 */
unsigned int Fakebeta_::update(const Trimat_& Dista, const Polymer_& Polymer)
{
    if (Dista.rno()!=Polymer.len()+2)
    {
	cerr<<"\n? Fakebeta_::update(): Dista:Polymer size mismatch\n";
	return(Lambda.len());	// unchanged size
    }
    
    register unsigned int Rno=Dista.rno();
    Distab.set_size(Rno); Distb.set_size(Rno); 	// resize
    Lambda.len(Rno); Dhj.len(Rno);
    
    /* NOTE: Rno held the actual size of Dista which is the number of
     * residues + 2 for the N/C -termini. Now adjust Rno so that is
     * is the number of residues
     */
    Rno-=2;
    
    register unsigned int i, j;
    register float Daf, Dbf, Dcf, Defgh, Dij;
    
    // calc the Lambda[] and Dhj[] values
    make_lambda(Dista, Polymer);
    
    /* alpha[i]:beta[j] distances: the underlying logic is the same as
     * in DRAGON 3.x but the cycles are rearranged for efficiency.
     * Note the special treatment of the N- and C-termini.
     */
    
    // N-terminus (i==0)
    for (j=1; j<=Rno; j++)
    {
	Defgh=0.5F*((float)Dista[j-1][0]+(float)Dista[j+1][0])-0.25F*(float)Dista[j+1][j-1];
	Distab[0][j]=get_dist((float)Dista[j][0], Defgh, Dhj[j], Lambda[j]);
    }
    
    // middle of the chain
    for (i=1; i<=Rno; i++)
    {
	Distab[i][i]=Polymer.abdist(i-1);	    // just copy prescribed value
	for (j=1; j<=Rno; j++)
	{
	    if (i==j) continue;
	    
	    Daf=float((i>=j-1)? Dista[i][j-1]: Dista[j-1][i]);
	    Dcf=float((i>=j+1)? Dista[i][j+1]: Dista[j+1][i]);
	    Defgh=0.5F*(Daf+Dcf)-0.25F*(float)Dista[j+1][j-1];
	    Dbf=float((i>=j)? Dista[i][j]: Dista[j][i]);
	    Distab[i][j]=get_dist(Dbf, Defgh, Dhj[j], Lambda[j]);
	}
    }
    
    // C-terminus (i==Rno+1 automatically here)
    for (j=1; j<=Rno; j++)
    {
	Defgh=0.5F*((float)Dista[i][j-1]+(float)Dista[i][j+1])-0.25F*(float)Dista[j+1][j-1];
	Distab[i][j]=get_dist((float)Dista[i][j], Defgh, Dhj[j], Lambda[j]);
    }
    
    // beta[i]:beta[j] distances
    for (i=2; i<=Rno; i++)
	for (j=1; j<i; j++)
	{
	    Dij=0.5F*((float)Distab[i-1][j]+(float)Distab[i+1][j])-0.25F*(float)Dista[i+1][i-1];
	    Distb[i][j]=get_dist((float)Distab[i][j], Dij, Dhj[i], Lambda[i]);
	}
    return(Rno);
}
// END of update()

/* beta_xyz(): generates the fake C-beta coordinates from the C-alpha coordinates
 * stored in Xyz and puts the result into Beta.
 */
void Fakebeta_::beta_xyz(const Points_& Xyz, const Polymer_& P, Points_& Beta)
{
    unsigned int Rno=P.len();	// all dims are checked beforehand
    Vector_ H(Xyz.dim());
    register unsigned int i;
    double Dbh, Dbj;
    
    for (i=1; i<=Rno; i++)
    {
        // get the midpoint coords and the B-H distance
	H=(Xyz[i-1]+Xyz[i+1]); H/=2.0;
        Dbh=diff_len2(Xyz[i], H);
	Dbj=P.abdist(i-1);    // prescribed distance
	
        if (Dbh==0.0 || Dbj==0.0)        /* i-1:i:i+1 angle is Pi; beta is on alpha */
            Beta[i]=Xyz[i];
        else
        {
            Dbh=sqrt(Dbj/Dbh);
	    Beta[i]=(Xyz[i]-H); Beta[i]*=Dbh;
            Beta[i]+=Xyz[i];
        }
    }
}
// END of beta_xyz()

// ---- Auxiliaries ----

/* make_lambda: given the C-alpha distances between A,B,C, and the
 * fake C-alpha/C-beta B-J distances in Seq[].Abdist, the H-J
 * distances are calculated for every C-beta [1..Rno], returned
 * in the vector Dhj (squared). Also the Lambda=
 * d(BJ)/d(BH) is calculated for each beta.
 * As there are no betas on the 0-th and Rno+1-th "alphas" (N/C termini), 
 * the corresponding entries in Dbj[], Dhj[] and Lambda[]
 * contain junk. Private
 */
void Fakebeta_::make_lambda(const Trimat_& Dista, const Polymer_& Polymer)
{
    register unsigned int i;
    register float ab, ac, bc, bh, bj, hj;
    unsigned int Rno=Polymer.len();   // dim mismatches were filtered before
    
    ab=Dista[1][0];
    Lambda[0]=Lambda[Rno+1]=1.0;    // no side chain on terminals 27/11/95
    
    for (i=1; i<=Rno; i++)
    {
	bc=Dista[i+1][i]; ac=Dista[i+1][i-1];
	bh=(ab+bc)/2.0F-ac/4.0F; 
	if (bh<0.0) bh=fabsf(bh);    // silent non-metric problem
	bj=sqrtf((float)Polymer.abdist(i-1)); hj=bj+sqrtf(bh); 
	Lambda[i]=(hj==0.0F || bj>hj)? 1.0F :bj/hj;
	Dhj[i]=hj*hj;
	ab=bc; 
    }
}
// END of make_lambda()

// ==== END OF METHODS Fakebeta.c++ ====
