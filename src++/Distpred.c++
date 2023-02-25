// ==== PROJECT DRAGON: METHODS Distpred.c++ ====

/* Interresidue distance prediction based on the conserved
 * hydrophobicity score.
 */

// SGI C++ 4.0, IRIX 5.3, 12. May 1995. Andras Aszodi

// ---- STANDARD HEADERS ----

#include <iostream.h>
#include <iomanip.h>

// ---- CLASS HEADER ----

#include "Distpred.h"

// ---- UTILITY HEADERS ----

#include "Pmest.h"

// ---- DEFINITIONS ----

#define DIST_BINNO 100
#define PARAMNO 3
#define MIN_DIST (0.0)
#define MAX_DIST (60.0)

// ---- STATIC INITIALISATION ----

Spl_ Distpred_::Idspl=Distpred_::init_idspl();

// ==== Distpred_ METHODS ====

// ---- Initialisation ----

/* init_idspl(): fills up the ideal C-alpha distance CDF spline with
 * the observed C-alpha distance distribution of a bunch of monomeric
 * proteins between 100 and 200 residues. See Aszodi & Taylor, J. Math. Chem.
 * for details.
 */
Spl_ Distpred_::init_idspl()
{
    Spl_ Ids(DIST_BINNO);
    
    /* In DRAGON 3.x, the ideal distribution was in a parameter
     * file. From 4.0 on, the values are hard-wired as they never
     * changed in 3.x.
     */
    Ids.x(0)=0.000000e+00; Ids.y(0)=0.000000e+00;
    Ids.x(1)=6.060606e-01; Ids.y(1)=0.000000e+00;
    Ids.x(2)=1.212121e+00; Ids.y(2)=0.000000e+00;
    Ids.x(3)=1.818182e+00; Ids.y(3)=0.000000e+00;
    Ids.x(4)=2.424242e+00; Ids.y(4)=0.000000e+00;
    Ids.x(5)=3.030303e+00; Ids.y(5)=0.000000e+00;
    Ids.x(6)=3.636364e+00; Ids.y(6)=2.291029e-05;
    Ids.x(7)=4.242424e+00; Ids.y(7)=3.751561e-04;
    Ids.x(8)=4.848485e+00; Ids.y(8)=2.557362e-03;
    Ids.x(9)=5.454545e+00; Ids.y(9)=1.244602e-02;
    Ids.x(10)=6.060606e+00; Ids.y(10)=2.441378e-02;
    Ids.x(11)=6.666667e+00; Ids.y(11)=3.585175e-02;
    Ids.x(12)=7.272727e+00; Ids.y(12)=4.471803e-02;
    Ids.x(13)=7.878788e+00; Ids.y(13)=5.147943e-02;
    Ids.x(14)=8.484848e+00; Ids.y(14)=6.139100e-02;
    Ids.x(15)=9.090909e+00; Ids.y(15)=7.635428e-02;
    Ids.x(16)=9.696970e+00; Ids.y(16)=9.321053e-02;
    Ids.x(17)=1.030303e+01; Ids.y(17)=1.149209e-01;
    Ids.x(18)=1.090909e+01; Ids.y(18)=1.366170e-01;
    Ids.x(19)=1.151515e+01; Ids.y(19)=1.565546e-01;
    Ids.x(20)=1.212121e+01; Ids.y(20)=1.781104e-01;
    Ids.x(21)=1.272727e+01; Ids.y(21)=2.016049e-01;
    Ids.x(22)=1.333333e+01; Ids.y(22)=2.262248e-01;
    Ids.x(23)=1.393939e+01; Ids.y(23)=2.526175e-01;
    Ids.x(24)=1.454545e+01; Ids.y(24)=2.804965e-01;
    Ids.x(25)=1.515152e+01; Ids.y(25)=3.083869e-01;
    Ids.x(26)=1.575758e+01; Ids.y(26)=3.363718e-01;
    Ids.x(27)=1.636364e+01; Ids.y(27)=3.646030e-01;
    Ids.x(28)=1.696970e+01; Ids.y(28)=3.930404e-01;
    Ids.x(29)=1.757576e+01; Ids.y(29)=4.209967e-01;
    Ids.x(30)=1.818182e+01; Ids.y(30)=4.490074e-01;
    Ids.x(31)=1.878788e+01; Ids.y(31)=4.777570e-01;
    Ids.x(32)=1.939394e+01; Ids.y(32)=5.066927e-01;
    Ids.x(33)=2.000000e+01; Ids.y(33)=5.357429e-01;
    Ids.x(34)=2.060606e+01; Ids.y(34)=5.638682e-01;
    Ids.x(35)=2.121212e+01; Ids.y(35)=5.912059e-01;
    Ids.x(36)=2.181818e+01; Ids.y(36)=6.183861e-01;
    Ids.x(37)=2.242424e+01; Ids.y(37)=6.445439e-01;
    Ids.x(38)=2.303030e+01; Ids.y(38)=6.698942e-01;
    Ids.x(39)=2.363636e+01; Ids.y(39)=6.946573e-01;
    Ids.x(40)=2.424242e+01; Ids.y(40)=7.183151e-01;
    Ids.x(41)=2.484848e+01; Ids.y(41)=7.403891e-01;
    Ids.x(42)=2.545455e+01; Ids.y(42)=7.615983e-01;
    Ids.x(43)=2.606061e+01; Ids.y(43)=7.812726e-01;
    Ids.x(44)=2.666667e+01; Ids.y(44)=8.003053e-01;
    Ids.x(45)=2.727273e+01; Ids.y(45)=8.181438e-01;
    Ids.x(46)=2.787879e+01; Ids.y(46)=8.341896e-01;
    Ids.x(47)=2.848485e+01; Ids.y(47)=8.490899e-01;
    Ids.x(48)=2.909091e+01; Ids.y(48)=8.635205e-01;
    Ids.x(49)=2.969697e+01; Ids.y(49)=8.766281e-01;
    Ids.x(50)=3.030303e+01; Ids.y(50)=8.887648e-01;
    Ids.x(51)=3.090909e+01; Ids.y(51)=8.996873e-01;
    Ids.x(52)=3.151515e+01; Ids.y(52)=9.091635e-01;
    Ids.x(53)=3.212121e+01; Ids.y(53)=9.182561e-01;
    Ids.x(54)=3.272727e+01; Ids.y(54)=9.263319e-01;
    Ids.x(55)=3.333333e+01; Ids.y(55)=9.335372e-01;
    Ids.x(56)=3.393939e+01; Ids.y(56)=9.400237e-01;
    Ids.x(57)=3.454545e+01; Ids.y(57)=9.458802e-01;
    Ids.x(58)=3.515152e+01; Ids.y(58)=9.509032e-01;
    Ids.x(59)=3.575758e+01; Ids.y(59)=9.555025e-01;
    Ids.x(60)=3.636364e+01; Ids.y(60)=9.595547e-01;
    Ids.x(61)=3.696970e+01; Ids.y(61)=9.635125e-01;
    Ids.x(62)=3.757576e+01; Ids.y(62)=9.667543e-01;
    Ids.x(63)=3.818182e+01; Ids.y(63)=9.695866e-01;
    Ids.x(64)=3.878788e+01; Ids.y(64)=9.722499e-01;
    Ids.x(65)=3.939394e+01; Ids.y(65)=9.746727e-01;
    Ids.x(66)=4.000000e+01; Ids.y(66)=9.766859e-01;
    Ids.x(67)=4.060606e+01; Ids.y(67)=9.785875e-01;
    Ids.x(68)=4.121212e+01; Ids.y(68)=9.803573e-01;
    Ids.x(69)=4.181818e+01; Ids.y(69)=9.820526e-01;
    Ids.x(70)=4.242424e+01; Ids.y(70)=9.835132e-01;
    Ids.x(71)=4.303030e+01; Ids.y(71)=9.848678e-01;
    Ids.x(72)=4.363636e+01; Ids.y(72)=9.860075e-01;
    Ids.x(73)=4.424242e+01; Ids.y(73)=9.872332e-01;
    Ids.x(74)=4.484848e+01; Ids.y(74)=9.883215e-01;
    Ids.x(75)=4.545455e+01; Ids.y(75)=9.893582e-01;
    Ids.x(76)=4.606061e+01; Ids.y(76)=9.902746e-01;
    Ids.x(77)=4.666667e+01; Ids.y(77)=9.911910e-01;
    Ids.x(78)=4.727273e+01; Ids.y(78)=9.920358e-01;
    Ids.x(79)=4.787879e+01; Ids.y(79)=9.927861e-01;
    Ids.x(80)=4.848485e+01; Ids.y(80)=9.934706e-01;
    Ids.x(81)=4.909091e+01; Ids.y(81)=9.941350e-01;
    Ids.x(82)=4.969697e+01; Ids.y(82)=9.946905e-01;
    Ids.x(83)=5.030303e+01; Ids.y(83)=9.952461e-01;
    Ids.x(84)=5.090909e+01; Ids.y(84)=9.958017e-01;
    Ids.x(85)=5.151515e+01; Ids.y(85)=9.963344e-01;
    Ids.x(86)=5.212121e+01; Ids.y(86)=9.967496e-01;
    Ids.x(87)=5.272727e+01; Ids.y(87)=9.971563e-01;
    Ids.x(88)=5.333333e+01; Ids.y(88)=9.975572e-01;
    Ids.x(89)=5.393939e+01; Ids.y(89)=9.979066e-01;
    Ids.x(90)=5.454545e+01; Ids.y(90)=9.982130e-01;
    Ids.x(91)=5.515152e+01; Ids.y(91)=9.984965e-01;
    Ids.x(92)=5.575758e+01; Ids.y(92)=9.987056e-01;
    Ids.x(93)=5.636364e+01; Ids.y(93)=9.989490e-01;
    Ids.x(94)=5.696970e+01; Ids.y(94)=9.991924e-01;
    Ids.x(95)=5.757576e+01; Ids.y(95)=9.994559e-01;
    Ids.x(96)=5.818182e+01; Ids.y(96)=9.996134e-01;
    Ids.x(97)=5.878788e+01; Ids.y(97)=9.997652e-01;
    Ids.x(98)=5.939394e+01; Ids.y(98)=9.998854e-01;
    Ids.x(99)=6.000000e+01; Ids.y(99)=1.000000e+00;
    // goes up to 60 A which is MAX_DIST
    
    Ids.fit_spl();
    return(Ids);
}
// END of init_idspl()

/* init_par(): inits the parameter vector to the parameter
 * values in the J. Math. Chem. paper.
 */
Vector_ Distpred_::init_par()
{
    Vector_ Par(PARAMNO);
    Par[0]=30.3; Par[1]=0.26; Par[2]=50.0;
    return(Par);
}
// END of init_par()

// ---- Transform parameter estimation ----

/* estim_params(): estimates the parameters of the hydrophobic score
 * transform function from the phobicity*conservation data in
 * Consphob. If Consphob is empty then no action is taken.
 */
void Distpred_::estim_params(const Array_<double>& Consphob)
{
    if (!Consphob.len())
    {
	cerr<<"\n? Distpred_::estim_params(): No sequences, no action\n";
	return;
    }
    
    Cdf_ Rawcdf=make_distr(Consphob);
    
    Vector_ W(DIST_BINNO), P(PARAMNO), Sd(PARAMNO);
    Trimat_ Correl(PARAMNO);
    float Q, Tstat95, Steplim=1e-6;
    int Itmax=200;
    
    W.set_values(1.0);	// uniform weight
    P=init_par();	// initial guess
    
    // do the estimation
    Q=nonlin11_reg(Rawcdf.x_vec(), Rawcdf.y_vec(), W, transform_hdist, 
		P, Sd, Correl, Tstat95,
		Itmax, Steplim, NLIN_TALK);
    cout<<"\nQ="<<Q<<", Stepno="<<Itmax<<", t-stat="<<Tstat95<<endl;
    cout.precision(3);
    cout<<"\nD=-"<<P[0]<<" * H ^"<<P[1]<<" + "<<P[2]<<endl<<endl;
    cout<<"Standard deviations:\n"<<Sd;
    cout<<"Correlation matrix:\n"<<Correl;
    
    Par=P;
}
// END of estim_params()

/* make_distr(): constructs and returns the cumulative distribution
 * function of the raw hydrophobic scores. Private
 */
Cdf_ Distpred_::make_distr(const Array_<double>& Consphob)
{
    // set up a vector for the cons*phob products
    unsigned int Rno=Consphob.len(); // already checked that >0 before this call
    double *Rawhd=new double [Rno]; // put cons*phob here
    double Hdmax=-1.0, Rhd;
    register unsigned int i, j;
    
    // fill up the raw score array, get maximal score
    for (i=0; i<Rno; i++)
    {
	Rhd=Consphob[i];
	if (Rhd<0.0) continue;	// skip negatives
	if (Rhd>Hdmax) Hdmax=Rhd;
	Rawhd[i]=Rhd;
    }

    Cdf_ Rawcdf(DIST_BINNO, 0.0, 2.0*fabs(Hdmax));	// set up CDF object
    for (i=0; i<Rno; i++)
    {
	if (Rawhd[i]<0.0) continue;	// skip negatives
	for (j=0; j<i; j++)
	{
	    if (Rawhd[j]<0.0) continue;
	    Rawcdf+=(Rawhd[i]+Rawhd[j]);	// add sum to CDF
	}
    }
    
    delete [] Rawhd;
    return(Rawcdf);
}
// END of make_distr()

/* transform_hdist(): this is the nonlinear function to be fitted.
 * H, a raw hydrophobic dist, is turned into a dist estimate
 * by D(H, P)=-P[0]*H^(P[1])+P[2],  Willie's modified empirical function.
 * F(H), the estimated CDF of the distribution of the hydrophobic
 * distances is:
 * F(H)=1-G(D(H)),  where G(D) is the observed CDF of the C-alpha
 * distances kept in the global spline Idspl. The transformation
 * takes into account that D'(H)<0 (mon.decr.). Private
 * Return value: F(H, P).
 */
double Distpred_::transform_hdist(double H, const Vector_& P)
{
    register double D, Gd;

    D=dist_phob(H, P);	// raw H into raw D
    
    // spline interpolation to find G(D(H))
    Gd=(D<MIN_DIST || D>MAX_DIST)? 0.0: Idspl.eval_spl(D);
    return(1.0-Gd);
}
// END of transform_hdist()

/* dist_phob(): returns the distance corresponding to the
 * raw hydrophobic estimate H or 0.0 if H<=0.0. It is assumed
 * that all P[i]>=0.0.
 */
double Distpred_::dist_phob(double H, const Vector_& P)
{
    if (H<0.0) return(0.0);
    if (H<DBL_EPSILON) return(P[2]);
    return(-P[0]*pow(H, P[1])+P[2]);
}
// END of dist_phob()

// ==== END OF METHODS Distpred.c++ ====
