// ==== PROJECT DRAGON: METHODS Iproj.c++ ====

/* The Hierarchic Inertial Projection. */

// SGI C++, IRIX 6.2, 14. Aug. 1996. Andris Aszodi 

// ---- STANDARD HEADERS ----

#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

// ---- MODULE HEADER ----

#include "Iproj.h"

// ---- UTILITY HEADERS ----

#include "Ql.h"
#include "Rsmdiag.h"
#include "Hirot.h"

// ---- TYPEDEFS AND PROTOTYPES ----

/* The following C-style constructs are necessary for a qsort() call. */
// struct for storing the distance quality for each cluster
typedef struct
{
    unsigned int ci, a0;	// cluster index,ctr pos.
    bool Flip;	// true if the cluster has been flipped
    int Detsign;    // sign of the rot matrix determinant
    float Q;	// dist quality
}
Rs_;

extern "C" int descend_rs(const void *P1, const void *P2);

// ==== METHODS ====

// ---- Constructors ----

/* Inits to perform projections on a point set made up of Resno points */
Iproj_::Iproj_(unsigned int Resno):
	Rno(Resno), Cluno(0), 
	Locals(Resno), Locdist(Resno), Clusters(NULL), 
	Ptclu(NULL), Ptoffs(NULL), Cluoffs(NULL), 
	Maxlocdim(0), Sksize(0), Diagshf(0.0)
{
    if (!Rno)
    {
	cerr<<"Iproj_("<<Rno<<", "<<Cluno<<"): Rno=1 assumed\n"; Rno=1;
	Locdist.set_size(Rno);
    }
    
    // set dimensions if necessary
    Locals.len_dim(Rno, Rno);
    make_clusters();
}

// ---- Size ----

/* set_size(): sets the calling object to Resno residues.
 * Works much in the same way as the corresponding constructor.
 * If Cl is 0 (the default), then the cluster number is calculated
 * within; any other number is interpreted as the required number
 * of clusters.
 */
void Iproj_::set_size(unsigned int Resno, unsigned int Cl)
{
    if (Rno==Resno && Cl==Cluno) return;	// no change
    
    Rno=Resno;
    if (!Rno)
    {
	cerr<<"Iproj_::set_size("<<Rno<<"): Rno=1 assumed\n"; Rno=1;
    }
    
    // zero cluster number indicates internally generated cluster nos.
    make_clusters(Cl);
    
    // set dimensions if necessary
    Locals.len_dim(Rno, Rno);
    Locdist.set_size(Rno);
    Maxlocdim=Sksize=0; Diagshf=0.0;
}
// END of set_size()

/* make_clusters(): constructs an array of bit-vectors (Cno all in all) which
 * store cluster membership information. All bit-vectors are Rno long, 
 * and the p-th bit is "true" in the ci-th bit-vector if the p-th point
 * is a member of the ci-th cluster. The union of the clusters is the
 * whole point set and their pairwise intersections are empty.
 * In the current scheme, "meshing" clusters are generated
 * along the chain, irrespective of the secstr layout, like this:-
 * [0] 100100100...
 * [1] 010010010...
 * [2] 001001001...
 * so that each cluster "covers" the whole chain. If Cno==0 (the default), 
 * Cluno is determined internally.
 * make_clusters(Clus): If an array of bit-vectors is given as an argument, 
 * then it is copied into Clusters and Cluno will be set to its length.
 * Return value: the number of clusters or 0 on error.
 */
unsigned int Iproj_::make_clusters(unsigned int Cno)
{
    if (Clusters!=NULL) delete [] Clusters;
    
    // determine cluster number
    if (!Cno)
    {
	Cluno=Rno/25+1;
	if (Cluno<=1) Cluno=2;
    }
    else Cluno=Cno;
    
    Clusters=new Bits_ [Cluno];
    if (Clusters==NULL)
    {
	cerr<<"\n! Iproj_::make_clusters(): Out of memory\n";
	Cluno=0;
	return(0);
    }
    
    // init each cluster
    unsigned int i, j;
    for (i=0; i<Cluno; i++)
    {
	Clusters[i].len(Rno); 
	Clusters[i].set_values(false);
    }
    for (j=0; j<Rno; j++)
    {
	i=j % Cluno;
	Clusters[i].set_bit(j);
    }
    Imoms.len_dim(Cluno, Rno);
    make_offsets();

    return(Cluno);
}

unsigned int Iproj_::make_clusters(const Array_<Bits_>& Clus)
{
    /* Various health checks: Clus must be non-empty,
     * all bitvectors must have the same length (Rno), their
     * union must be fully activated, pairwise intersection
     * must be empty. If any of the checks fails, the clusters
     * will be generated internally by make_clusters() above.
     * The 3 Bits_ objects are init'd at the beginning because
     * SGI C++ 6.2 didn't like the bypassed initialisations
     */
    Bits_ Union(Rno), Overlap(Rno), Smalls(Rno);
    
    if (!Clus.len())
    {
	cerr<<"\n? Iproj_::make_clusters(Clus): No clusters";
	goto BADCLUS;
    }
    
    register unsigned int i, j;
    Union.set_values(false);
    for (i=0; i<Clus.len(); i++)
    {
	if (Clus[i].len()!=Rno)	// length mismatch
	{
	    cerr<<"\n? Iproj_::make_clusters(Clus): Cluster "<<i
		<<" has length "<<Clus[i].len()<<" instead of "<<Rno;
	    goto BADCLUS;
	}
	Union|=Clus[i];
	
	// check pairwise overlaps
	for (j=0; j<i; j++)
	{
	    Overlap=Clus[i] & Clus[j];
	    if (Overlap.on_no())
	    {
		cerr<<"\n? Iproj_::make_clusters(Clus): Cluster "<<i
		    <<" overlaps with cluster "<<j;
		goto BADCLUS;
	    }
	}	// for j
    }	    // for i
    
    if (Union.on_no()!=Rno) // no full coverage
    {
	cerr<<"\n? Iproj_::make_clusters(Clus): Clusters don't cover full point set";
	goto BADCLUS;
    }
    
    /* End of health checks */
    
    // merge small clusters into one
    register unsigned int Smallno;
    for (i=Smallno=0; i<Clus.len(); i++)
	if (Clus[i].on_no()<=3)
	{
	    Smalls|=Clus[i];
	    Smallno++;
	}
    
    // will be separate cluster if contains at least 6 points
    Cluno=(Smalls.on_no()>5)? Clus.len()-Smallno+1: Clus.len()-Smallno;

    // store external clusters
    if (Clusters!=NULL) delete [] Clusters;
    Clusters=new Bits_ [Cluno];
    if (Clusters==NULL)
    {
	cerr<<"\n! Iproj_::make_clusters(Clus): Out of memory\n";
	Cluno=0;
	return(0);
    }
    
    for (i=j=0; i<Clus.len(); i++)
    {
	if (Clus[i].on_no()<=3) continue;   // don't save smalls
	Clusters[j++]=Clus[i];
    }
    if (Smalls.on_no()>5) Clusters[j]=Smalls;	// last cluster is the joint small
    else Clusters[0]|=Smalls;	// merge smalls with first

    Imoms.len_dim(Cluno, Rno);
    make_offsets();
    return(Cluno);
    
    /* End of normal execution. If any of the cluster checks failed,
     * an ugly "goto" brings us here where the default clusters
     * are constructed (and also the error messages are finished).
     */
    BADCLUS:
    make_clusters();
    cerr<<" ("<<Cluno<<" generated internally)\n";
    return(Cluno);

}
// END of make_clusters()

/* make_offsets(): constructs and fills up 3 uint arrays which
 * contain various offsets. Ptclu[i] is the no. of the cluster
 * containing the i:th point. Ptoffs[i] returns the column index
 * in Abprods[][] (cf. ctr_prod()) corresponding to the i:th point.
 * Cluoffs[ci] returns the column index in Abprods for the ci:th
 * cluster. These index arrays are meant to speed up the ctr_prod()
 * method. Private
 */
void Iproj_::make_offsets()
{
    /* resize the arrays if they were allocated */
    if (Ptclu!=NULL) delete [] Ptclu;
    Ptclu=new unsigned int [Rno];
    
    if (Ptoffs!=NULL) delete [] Ptoffs;
    Ptoffs=new unsigned int [Rno];

    if (Cluoffs!=NULL) delete [] Cluoffs;
    Cluoffs=new unsigned int [Cluno];

    register unsigned int i, ci, k;
    
    // fill up Ptclu
    for (i=0; i<Rno; i++)
    {
	// find the cluster containing the i:th point
	for (ci=0; ci<Cluno && !Clusters[ci].get_bit(i); ci++);
	Ptclu[i]=ci;
    }
    
    // init the Abprods offset arrays
    for (ci=k=0; ci<Cluno; ci++)
    {
	Cluoffs[ci]=k++;
	for (i=0; i<Rno; i++)
	    if (Clusters[ci].get_bit(i))
		Ptoffs[i]=k++;
    }
}
// END of make_offsets()

// ---- Projections ----

/* full_project(): performs the Hierarchic Inertial Projection
 * on a point set. Dist holds the squared interpoint distances.
 * The projections will use an Evfract-th fraction of
 * the sum of all positive eigenvalues and will project into a
 * Dim<Oldim-dimensional Euclidean space. The Cartesian coordinates
 * will be put into Xyz (the sizes of which is set within).
 * Return value: the new dimension.
 */
unsigned int Iproj_::full_project(Trimat_& Dist,
	double Evfract, unsigned int Oldim, Points_& Xyz)
{
    // check sizes
    if (Dist.rno()!=Rno)
    {
	cerr<<"\n? Iproj_::full_project(): Size mismatch\n";
	return(0);
    }
    
    Trimat_ Metric(Rno);
    trineq_filter(Dist, Metric);	// full tri balancing

    // don't do local projections if there's just one cluster
    if (Cluno==1)
    {
	Xyz.mask(true);	    // switch everybody ON
	return(metric_project(Metric, Evfract, 3, Oldim, Xyz));
    }
    
    // "crushed" local embedding of clusters
    cluster_project(Dist, Oldim);

    // skeleton embedding
    unsigned int Df=skel_project(Dist, Evfract, Oldim, Xyz);

    return(Df); // full proj dim returned
}
// END of full_project()

/* cluster_project(): given a distance matrix Dist,
 * the local cluster coordinates and moments of
 * inertia are updated. Also sets the maximal local dimension
 * and the skeleton size.
 * Dist is modified: after local projection, the new local distances
 * (which are metric now) are put back into it.
 * No action is taken when invoked on 1-cluster sets.
 */
void Iproj_::cluster_project(Trimat_& Dist, unsigned int Oldim)
{
    if (Cluno<=1) return;   // don't process single clusters
    
    Trimat_ Dloc, Mloc;   // local dist and metric matrix
    register unsigned int ci, D, Embed;
    
    // local embedding of clusters
    Sksize=Maxlocdim=0;	    // reset skeleton size and max. local dimension
    Diagshf=0.0;	// erase any memories of prev triangle filterings
    Locdist.set_values();   // zero the local Euclidean distance matrix
    for (ci=0; ci<Cluno; ci++)
    {
	// 1-point clusters (zero local dim) are special
	if (Clusters[ci].on_no()==1)
	{
	    Locals.mask(Clusters[ci]);   // mask it
	    Locals[0].set_values();	    // local center is the null-vector
	    continue;
	}
	
	sub_matrix(Dist, Clusters[ci], Dloc);
	dist_metric(Dloc, Mloc);
	Locals.mask(Clusters[ci]);	// select the ci-th cluster
	Embed=Mloc.rno();
	if (Oldim<Embed) Embed=Oldim;
	
	D=metric_project(Mloc, 1.0, 1, Embed, Locals, &(Imoms[ci]));
	
	make_locdist();	// create local Euclidean distances
	Sksize+=D;    // sum dimensions for skeleton size
	if (D>Maxlocdim) Maxlocdim=D;	// get max. local dimension
    }
    apply_locdist(Dist);
    Sksize+=Cluno;  // the final skeleton size
}
// END of cluster_project()

/* skel_project(): Performs the skeleton projection. The "skeleton" is made up
 * of the centroids of the clusters plus "inertial points" sitting on
 * the inertial axes of the clusters, a moment of inertia away from
 * the local centroid. The method assumes that the local coordinates
 * of the clusters have been determined beforehand by cluster_project().
 * Dm contains a distance matrix which
 * has modified inter-cluster entries but the intra-cluster distances
 * correspond to the previous local projections done in full_project().
 * The modified dist matrix is then first converted into a metric
 * matrix and then the skeleton projection is performed. Xyz will
 * contain the full Euclidean embedding on return.
 * Return value: the embedding dimension.
 */
unsigned int Iproj_::skel_project(Trimat_& Dm,
	double Evfract, unsigned int Oldim, Points_& Xyz)
{
    if (Cluno<=1) return(Oldim);    // don't do anything in 1-cluster case
    
    static Trimat_ Metric;
    dist_metric(Dm, Metric);
    
    // fill up the skeleton metric matrix
    Skmet.set_size(Sksize);

    make_skmet(Metric);

    // embed the skeleton
    Skxyz.len_dim(Sksize, Oldim);
    unsigned int Dim;

    // disallow flat skeletons
    Dim=metric_project(Skmet, Evfract, 
	(Maxlocdim>=3)? Maxlocdim: 3, Oldim, Skxyz);

    // put the local flesh onto the skeleton
    flesh_skel(Dm, Xyz);
    
    return(Dim);
    
}
// END of skel_project()

// ---- Reconstruction ----

/* flesh_skel(): puts the local structures in Locals onto the
 * Euclidean skeleton Skxyz so that the final structure is
 * returned in Xyz. Private
 */
void Iproj_::flesh_skel(const Trimat_& Dist, Points_& Xyz)
{
    register unsigned int ci, wi, a0, Da, p, i, Dim=Skxyz.dim(), Fno;
    Matrix_ R;
    Vector_ Iv(Dim);
        
    static Hirot_ Hr;
    static Array_<Points_> Distorts, Ideals;	// for each clu
    static Points_ W;	// weights
    Rs_ *Rss=new Rs_ [Cluno];	// RMS values for the cluster fits
    
    Xyz.len_dim(Rno, Dim);	// activate all Xyz
    Distorts.len(Cluno); Ideals.len(Cluno);
    W.len(Cluno); W.mask(true);
    
    const double HALFVAL=0.1;
    double Lendiff;
    int Detsign;

    Xyz.len_dim(Rno, Dim);	// activate all Xyz
    for (ci=wi=a0=0; ci<Cluno; ci++)
    {
	// check the cluster size: don't attempt rotation if it is 1
	Xyz.mask(Clusters[ci]);
	if (Xyz.active_len()==1)    // 1-point cluster
	{
	    Xyz[0]=Skxyz[a0];	// the whole thing is just a centroid
	    continue;
	}
	
	// make the local/global "rotation" matrix
	Locals.mask(Clusters[ci]);	// activate ci-th cluster only
	
	Da=Locals.dim();
	R.set_size(Dim, Da);
	
	Distorts[ci].len_dim(Da, Dim);
	Ideals[ci].len_dim(Da, Dim);
	W[ci].dim(Da);
	
	// make up ideal and distorted frames
	for (p=0; p<Da; p++)
	{
	    Distorts[ci][p]=Skxyz[a0+p+1]-Skxyz[a0];   // local p-th axis
	    Ideals[ci][p].set_values();
	    Lendiff=Ideals[ci][p][p]=Imoms[ci][p];  // p-th coord is the ideal length
	    Lendiff=fabs(Distorts[ci][p].vec_len()-Lendiff)/Lendiff;
	    W[ci][p]=(HALFVAL/(HALFVAL+Lendiff));	    // weight is high if rel. length diff is small
	}
	
	// rotate ideal frame onto distorted (with flip: improper rotation)
	Hr.best_rotflip(Ideals[ci], Distorts[ci], W[ci]);
	
	// get the sign of the determinant of the rotation matrix
	if (Detsign=Hr.det_sign())	// assignment intentional
	{	// save index etc. for this cluster
	    Rss[wi].ci=ci; Rss[wi].a0=a0;
	    Rss[wi].Flip=false; Rss[wi].Detsign=Detsign;
	    ++wi;
	}
	
	for (p=0; p<Da; p++)
	{
	    Iv=Hr.rot_matrix()*Ideals[ci][p];	    // new local p-th axis
	    Iv.vec_norm();
	    R.col(Iv, p);
	}
	
	// apply the local-->global transform to points in local cluster
	for (i=0; i<Locals.active_len(); i++)
	    Xyz[i]=R*Locals[i]+Skxyz[a0];

	a0+=Da+1;    
    }
    Xyz.mask(true);
    
    Fno=wi;	// the number of flippable ("non-flat") clusters
    if (!Fno)	// don't bother w/ corrections if no flips can be done...
    {
	delete [] Rss;
	return;
    }
    
    // see if flipping the clusters improves the overall structure
    float Q, Qflip;
    Points_ Xyzflip=Xyz;    // temp copy for flipping local clus
    register unsigned int Flipno;
    
    /* NOTE: redoing the rotations is terribly inefficient:
     * if it works, then think about cleverly arranging the -1.0
     * multiplications to save time!!!
     */
    do
    {
	/* Check the quality of each cluster by calculating the
	 * dist difference between its points and all other points
	 * in the other clusters. Ideally, this should be 0.0
	 */
	for (wi=0; wi<Fno; wi++)
	{
	    if (Rss[wi].Flip) continue;	// don't bother, was flipped already
	    // get the cluster quality
	    Rss[wi].Q=clu_qual(Clusters[Rss[wi].ci], Xyzflip, Dist);
	}
	qsort(Rss, Fno, sizeof(Rs_), descend_rs);	// sort in descending order

	for (Flipno=wi=0; wi<Fno; wi++)
	{
	    // start with the highest RMS unflipped
	    if (Rss[wi].Flip) continue;
	    
	    ci=Rss[wi].ci; a0=Rss[wi].a0;
	    Q=Rss[wi].Q;
	    
	    Locals.mask(Clusters[ci]);	// activate ci-th cluster only
	    Da=Locals.dim();
	    R.set_size(Dim, Da);
	    
	    // invert the shortest axis if the orig. transform was "pure"
	    if (Rss[wi].Detsign>0)
		Ideals[ci][Da-1][Da-1]*=(-1.0);
	    
	    // "pure" rotation onto distorted
	    Hr.best_rot(Ideals[ci], Distorts[ci], W[ci]);
	    
	    for (p=0; p<Da; p++)
	    {
		Iv=Hr.rot_matrix()*Ideals[ci][p];	    // new local p-th axis
		Iv.vec_norm();
		R.col(Iv, p);
	    }
	    
	    // apply the flip transform to the copy
	    Xyzflip.mask(Clusters[ci]);
	    for (i=0; i<Locals.active_len(); i++)
		Xyzflip[i]=R*Locals[i]+Skxyz[a0];
	    
	    // check if the flipping improved matters
	    Xyzflip.mask(true);
	    Qflip=clu_qual(Clusters[ci], Xyzflip, Dist);
	    if (Qflip<Q)
	    {
		Flipno++;	// got better
		Rss[wi].Q=Qflip; Rss[wi].Flip=true;
		break;
	    }
	    else	// put back original vectors
	    {
		Xyz.mask(Clusters[ci]);
		Xyzflip.mask(Clusters[ci]);
		for (p=0; p<Clusters[ci].on_no(); p++)
		    Xyzflip[p]=Xyz[p];
		Xyzflip.mask(true);
		if (Rss[wi].Detsign>0)
		    Ideals[ci][Da-1][Da-1]*=(-1.0);	// flip back if necessary
	    }
	}
	Xyz.mask(true);
	if (Flipno) Xyz=Xyzflip;	// was modified, put back
    }
    while (Flipno);
    delete [] Rss;
}
// END of flesh_skel()

/* clu_qual(): given a set of points Xyz (assumed to be fully active)
 * and a target distance matrix Dist, check how different the distances
 * between the members of the cluster Clu are. Returns a quality
 * value (0.0 for perfect agreement). The sizes of Xyz and Dist
 * are assumed to be compatible. Static private
 */
float Iproj_::clu_qual(const Bits_& Clu, const Points_& Xyz, const Trimat_& Dist)
{
    register float Q=0.0, D;
    register unsigned int i, j, Pno=0, Size=Dist.rno();
    
    for (i=0; i<Size; i++)
	for (j=0; j<i; j++)
	{
	    if (Clu.get_bit(i)!=Clu.get_bit(j))
	    {
		D=diff_len2(Xyz[i], Xyz[j]);
		D-=Dist[i][j];
		Q+=fabs(D);
		Pno++;
	    }
	}
    return(Pno? Q/Pno: 0.0);
}
// END of clu_qual()

/* descend_rs(): aux function to sort the Rs_ records in descending
 * Relchg order. This function has C linkage.
 */
int descend_rs(const void *P1, const void *P2)
{
    const Rs_ *Wp1=(const Rs_*)P1, *Wp2=(const Rs_*)P2;
    if (Wp2->Q>Wp1->Q) return(1);
    else if (Wp2->Q<Wp1->Q) return(-1);
    else return(0);
}
// END of descend_rs()

// ---- Private projections ----

/* metric_project(): projects the metric matrix Metric into Euclidean space
 * with less than Oldim dimensions. Only the Evfract-th fraction of the
 * positive eigenvalues are used. However, the embedding dimension will not
 * be less than Mindim. 0<Mindim<=Oldim<Rno always holds true.
 * If Oldim<Mindim, then Mindim is set to Oldim.
 * Xyz is supposed to have been appropriately masked; the point dimensions
 * will be adjusted within. 
 * Return value: the new dimension or 0 if Xyz doesn't have enough 
 * active points. Also returns the sqrt of moments of inertia in *Moms
 * if Moms!=NULL (NULL is the default). Private
 */
unsigned int Iproj_::metric_project(const Trimat_& Metric,  
	double Evfract, unsigned int Mindim, unsigned int Oldim, 
	Points_& Xyz, Vector_ *Moms)
{
    unsigned int i, j, Size=Metric.rno();

    // diagonalise metric matrix
    static Vector_ Eval;
    static Sqmat_ Evec;
    
    if (Size>Xyz.active_len())	// paranoia
    {
	cerr<<"? metric_project(): free point no. "
	    <<Xyz.active_len()<<"<"<<Size<<", cannot project\n";
	return(0);  // *Moms not modified
    }
    Eval.dim(Size);
    Evec.set_size(Size);
    
    // modify dimension limits if necessary
    if (Oldim>Size) Oldim=Size;
    if (!Mindim) Mindim=1;
    if (Mindim>=Oldim) Mindim=Oldim-1;
    
    /* Diagonalisation: if the old dimension is less than fourth
     * of Size, then the Rsm object is used which generates all
     * eigenvalues and a chosen set of eigenvectors. Otherwise, 
     * all eigenvectors are made: this is faster, say the books.
     */
    Rsmdiag_ Rsm;
    bool Somevec=bool(4*Oldim<=Size);
    
    if (Somevec) Rsm.get_evals(Metric, Eval); // all eigenvalues
    else eigen_ql(Metric, Eval, Evec);	// total diagonalisation
    
    // subtract Diagshf from the eigenvalues
    for (i=0; i<Oldim; i++)
	Eval[i]-=Diagshf;
    
    // find the new dimension
    unsigned int Firstpos, Dim;
    
    // get smallest positive eigenvalue [1..N] style
    for (Firstpos=Oldim; Firstpos>0 && Eval[Firstpos-1]<FLT_EPSILON; Firstpos--);
    if (Firstpos<Mindim)
    {
	cerr<<"? metric_project(): only "<<Firstpos<<" positive eigenvalue";
	if (Firstpos!=1) cerr<<'s';
	cerr<<", cannot embed in Mindim="<<Mindim<<endl;
	Mindim=Firstpos;
    }
    
    // make up cumulative sums of eigenvalues
    Vector_ Sumeval(Firstpos);
    Sumeval[0]=Eval[0];
    for (i=1; i<Firstpos; i++)
	Sumeval[i]=Eval[i]+Sumeval[i-1];

    // determine no. of dimensions for Equ-th of total
    double Qu=Sumeval[Firstpos-1]*Evfract;
    for (Dim=1; Dim<=Firstpos && Sumeval[Dim-1]<Qu; Dim++);

    if (Dim>=Oldim) Dim=Oldim-1;	// force embedding in lower dimensions
    if (Dim<Mindim) Dim=Mindim;		// embed into at least Mindim dimensions

    // get the first Dim eigenvectors
    if (Somevec) Rsm.get_evecs(Dim, Evec);
    
    // generate Cartesian coordinates
    Xyz.dim(Dim); 
    for (j=0; j<Dim; j++)
    {
	Eval[j]=sqrt(Eval[j]);	// square root of largest evals
	if (Moms!=NULL) (*Moms)[j]=Eval[j];   // save sqrts moments of inertia
	
	// flip eigenvectors so that Evec[0][j]>0.0
	if (Evec[0][j]<0.0)
	    for (i=0; i<Size; i++) Evec[i][j]*=(-1.0);
    }
    
    for (i=0; i<Size; i++)
    {
	for (j=0; j<Dim; j++)
	    Xyz[i][j]=Eval[j]*Evec[i][j];	// column eigenvectors
    }
    
    return(Dim);
}
// END of metric_project()

// ---- Scalar products ----

/* make_skmet(): constructs the skeleton metric matrix Skmet, i.e. a matrix
 * that holds the scalar products of the cluster centroid vectors
 * and the inertial points. Skmet is assumed to have been allocated
 * to the correct size prior to the call.
 * Metric is the overall metric matrix of the individual points, 
 * Locals and Imoms hold the local inertial coordinates and moments of
 * inertia for all points. Locals is a "maskable array".
 * Imoms is supposed to be in an all-active
 * state (not checked). Private
 */
void Iproj_::make_skmet(const Trimat_& Metric)
{
    static Matrix_ Abprods;	// centroid scalar products
    static Array_<double> Momscal;	// moment scaling
    
    Abprods.set_size(Cluno, Rno+Cluno);
    Momscal.len(Sksize);
    ctr_prod(Metric, Abprods);	// get centroid scalar products

    register unsigned int ci, cj, ix, a0, b0, Da, Db, Pa, Pb, Na, Nb, p, q, ap, aq, bq;
    register double AA, AB, Spa0;
    Matrix_ Loca, Aibj, Sptq;	// temp matrices for mixed inertial scalprods
    Matrix_ *Locds=new Matrix_ [Cluno];    // local coords extracted for speed
    
    for (ci=a0=Pa=0; ci<Cluno; ci++)
    {
	// INTRA-cluster
	AA=Skmet[a0][a0]=Abprods[ci][Pa];   // <a0|a0>
	Locals.mask(Clusters[ci]);
	Na=Clusters[ci].on_no();	// dim and no. of pts in [ci]
	Da=(Na==1)? 0: Locals.dim();	// 1-member clusters are 0-dimensional
	
	for (p=0, ap=a0+1; p<Da; p++, ap++)
	{
	    // <s'p|a0>: norm s'p to sqrt of inert. mom.
	    Momscal[ap]=1.0/Imoms[ci][p];
	    Spa0=Skmet[ap][a0]=iv_ctrprod(Abprods, ci, Pa, p)*Momscal[ap];
	    
	    // <sp|sp>
	    Skmet[ap][ap]=Imoms[ci][p]*Imoms[ci][p]+2.0*Spa0+AA;
	    // <sp|sq>
	    for (aq=a0+1; aq<ap; aq++)
		Skmet[ap][aq]=Spa0+Skmet[aq][a0];   // AA has been added
	    
	    // <sp|a0>
	    Skmet[ap][a0]+=AA;	// AA added here before the [ap][aq] ^ there
	}

	// extract the local coords of [ci]-th cluster if not 0-dim
	if (Da)
	{
	    Locds[ci].set_size(Na, Da);
	    for (ix=0; ix<Na; ix++)
		Locds[ci].row(Locals[ix], ix);
	    Loca=Locds[ci].get_transpose();	    // transpose needed for ci/cj stuff
	}
	
	// INTER-cluster (between [ci] and [cj])
	for (cj=b0=Pb=0; cj<ci; cj++)
	{
	    AB=Skmet[a0][b0]=Abprods[ci][Pb];	// <a0|b0>
	    Locals.mask(Clusters[cj]);
	    Nb=Clusters[cj].on_no(); // dim and no. of pts in [cj]
	    Db=(Nb==1)? 0: Locals.dim(); // 1-point clusters are 0-dimensional
	    
	    // fill up the <a0|tq> row with <a0|t'q>: uses cj-th mask
	    for (q=0, bq=b0+1; q<Db; q++, bq++)
		Skmet[a0][bq]=iv_ctrprod(Abprods, ci, Pb, q)*Momscal[bq];
	    
	    // fill up the <sp|b0> column with <s'p|b0>: needs ci-th mask
	    Locals.mask(Clusters[ci]);
	    for (p=0, ap=a0+1; p<Da; p++, ap++)
		Skmet[ap][b0]=iv_ctrprod(Abprods, cj, Pa, p)*Momscal[ap];

	    // make the <sp|tq> rectangular region
	    Locals.mask(true);	// access will be via get_bit() from now on
	    Aibj.set_size(Na, Nb);	// <a'i|b'j>
	    
	    // prepare for <s'p|t'q> products if Da and Db both >0
	    if (Da>0 && Db>0)
	    {
		aibj_prod(Metric, Abprods, ci, cj, Pa, Pb, Aibj);
		Sptq=Loca*(Aibj*Locds[cj]);	    // this is <s'p|t'q>,not norm
	    }

	    // now calc global <sp|tq> and adjust the <sp|b0> column
	    for (p=0, ap=a0+1; p<Da; p++, ap++)
	    {
		for (q=0, bq=b0+1; q<Db; q++, bq++)
		    Skmet[ap][bq]=Sptq[p][q]*Momscal[ap]*Momscal[bq]    // norm now
			    +Skmet[ap][b0]+Skmet[a0][bq]+AB;
		Skmet[ap][b0]+=AB;
	    }
	    
	    // adjust <a0|tq> row
	    for (bq=b0+1; bq<Db+b0+1; bq++) Skmet[a0][bq]+=AB;
	    
	    // next step
	    b0+=Db+1; Pb+=Nb+1;
	}	// for cj
	
	// next step
	a0+=Da+1; Pa+=Na+1;
    }	    // for ci
    
    delete [] Locds;
}
// END of make_skmet()

/* ctr_prod(): constructs a rectangular matrix that holds the scalar products
 * of the cluster centroid vectors and the individual point vectors.
 * There are Cluno rows and Rno+Cluno columns and the layout is:-
 * <a0|a0> <a0|a1> ... <a0|an> | <a0|b0> <a0|b1> ... <a0|bm> | ....
 * <b0|a0> <b0|a1> ... <b0|an> | <b0|b0> <b0|b1> ... <b0|bm> | ....
 * ...
 * Input: the overall metric matrix Metric.
 * Output: the Abprods matrix (size adjusted silently if necessary). Private
 */
void Iproj_::ctr_prod(const Trimat_& Metric, Matrix_& Abprods) const
{
    register unsigned int ci, cj, i, j, ki, kj, kci, kcj, Colno=Rno+Cluno;
    register double Temp;
    
    Abprods.set_size(Cluno, Colno);   // adjust size if necessary
    Abprods.set_values();   // zero
    
    /* Fill up the Abprods matrix. Scan the metric matrix and 
     * sum the entries in the appropriate places. The three index
     * arrays (Ptclu, Ptoffs, Cluoffs) will provide the necessary
     * column indices. See make_offsets() 
     */
    for (i=0; i<Rno; i++)
    {
	ci=Ptclu[i];	// i is in the ci:th cluster
	ki=Ptoffs[i];	// and the <.0|ai> sum is in the ki-th column
	kci=Cluoffs[ci];
	
	Temp=Metric[i][i];
	Abprods[ci][kci]+=Temp;	    // <a0|a0>
	Abprods[ci][ki]+=Temp;	// <a0|ai>

	for (j=0; j<i; j++)
	{
	    cj=Ptclu[j];
	    kj=Ptoffs[j];
	    kcj=Cluoffs[cj];
	    
	    Temp=Metric[i][j];
	    
	    Abprods[ci][kcj]+=Temp;	// <a0|b0>
	    Abprods[cj][kci]+=Temp;	// <b0|a0> inefficient...
	    Abprods[ci][kj]+=Temp;	// <a0|bj>
	    Abprods[cj][ki]+=Temp;	// <b0|ai>
	}
    }

    // normalisation
    for (ci=0; ci<Cluno; ci++)
    {
	if (Clusters[ci].on_no()<=1) continue;
	
	Temp=1.0/Clusters[ci].on_no();	// membership
	for (kj=0; kj<Colno; kj++)    // divide the ci:th row
	    Abprods[ci][kj]*=Temp;
	kci=Cluoffs[ci];
	for (cj=0; cj<Cluno; cj++)
	    Abprods[cj][kci]*=Temp; // divide ci's column (kci:th)
    }
}
// END of ctr_prod()

/* iv_ctrprod(): returns the scalar product of the centroid of the
 * Ctridx-th cluster and the Coord-th local axis of inertia of the current
 * cluster. Abprods contains the scalar products of the position vectors
 * and the centroids (for the layout, see the comments of ctr_prod() ), 
 * Locals holds the local cluster coordinates, it is assumed to have been
 * masked appropriately before the call.
 * Coffs marks the column of Abprods where the Cluidx-th
 * cluster starts (with the centroid). Private
 */
double Iproj_::iv_ctrprod(const Matrix_& Abprods,
	unsigned int Ctridx, unsigned int Coffs, unsigned int Coord) const
{
    register unsigned int ic;
    register double A0b0=Abprods[Ctridx][Coffs], Temp=0.0;

    for (ic=0; ic<Locals.active_len(); ic++)
	Temp+=Locals[ic][Coord]*(Abprods[Ctridx][Coffs+ic+1]-A0b0);
    return(Temp);
}
// END of iv_ctrprod()

/* aibj_prod(): constructs the matrix Aibj (size is set before the call)
 * the [i][j]-th element of which is <ai|bj>-<a0|bj>-<ai|b0>+<a0|b0>.
 * <ai|bj> comes from Metric (the overall metric matrix), the rest
 * from Abprods. cluster A is the Aidx-th, 
 * cluster B is the Bidx-th. Aoffs and Boffs are the col idx offsets
 * for the centroids of A and B in Abprods, respectively. Private
 */
void Iproj_::aibj_prod(const Trimat_& Metric, const Matrix_& Abprods, 
	unsigned int Aidx, unsigned int Bidx, 
	unsigned int Aoffs, unsigned int Boffs, 
	Matrix_& Aibj) const
{
    register unsigned int i, ic, j, jc;
    register double Aib0, A0b0=Abprods[Aidx][Boffs];
    
    for (i=ic=0; i<Rno; i++)
    {
	if (!Clusters[Aidx].get_bit(i)) continue;
	
	Aib0=Abprods[Bidx][Aoffs+ic+1]-A0b0;	// <ai|b0>-<a0|b0>
	
	for (j=jc=0; j<Rno; j++)
	{
	    if (!Clusters[Bidx].get_bit(j)) continue;
	    
	    Aibj[ic][jc]=Metric(i, j)-Abprods[Aidx][Boffs+jc+1]-Aib0;
	    jc++;
	}
	ic++;
    }
}
// END of aibj_prod()

// ---- Triangle inequality balancing ----

/* trineq_filter(): smoothes triangle inequality violations in the
 * distance matrix Dist. Produces the corresponding metric matrix Metric.
 * This is an iterative routine that decides on the no. of iterations
 * internally if Tsmcyc=0 (the default). If Mass!=NULL, then 
 * Mass[i] will be interpreted as the mass of the i:th point
 * (for skel filters). Mass==NULL by default.
 * Return value: the actual no. of iterations done. Metric's dim is adjusted. Private
 */
unsigned int Iproj_::trineq_filter(Trimat_& Dist, Trimat_& Metric,
    unsigned int Tsmcyc, const double *Mass)
{
    const unsigned int TSM_FRAC=10, MIN_TSMCYC=1, MAX_TSMCYC=100;
    
    if (!Tsmcyc) Tsmcyc=Dist.rno()/TSM_FRAC;	// decided here
    if (Tsmcyc<MIN_TSMCYC) Tsmcyc=MIN_TSMCYC;
    if (Tsmcyc>MAX_TSMCYC) Tsmcyc=MAX_TSMCYC;
    
    unsigned int Itno, Tviol, Cviol;
    static Vector_ Cdist2;
    
    // construct Metric and smooth it until no violations are found
    Metric.set_size(Dist.rno());
    Cdist2.dim(Dist.rno());
    Diagshf=0.0;
    for (Itno=0; Itno<Tsmcyc; Itno++)
    {
	if (Mass==NULL)    // get distances from centroid
	    Cviol=centre_dist(Dist, Cdist2);	// all points are of unit mass
	else
	    Cviol=centre_dist(Dist, Mass, Cdist2);  // with different masses
	
	dist_metric(Dist, Cdist2, Metric);  // make metric matrix
	Tviol=trieq_bal(Metric);   // balance inequalities
	if (!Tviol && !Cviol) break;	    // success
	
	metric_dist(Metric, Dist);	    // get new dist matrix
    }
    return(Itno);
}
// END of trineq_filter()

/* trieq_bal: balances the triangle inequalities in Metric.
 * If its diagonal has non-negative values, then it is "shifted": 
 * twice the abs value of the worst non-negative diagonal element
 * is added to all diag elements. The shift value is stored in
 * Diagshf (for the projection). Then every off-diag entry is
 * checked for triangle inequality violations: being a scalar product,
 * cos(phi)=Metric[i][j]/sqrt(Metric[i][i]*Metric[j][j]), which should
 * fall in the range -1.0 ... +1.0. If not, then Metric[i][j] is scaled
 * a bit. 
 * Return value: the number of violations. Private
 */
unsigned int Iproj_::trieq_bal(Trimat_& Metric)
{
    const double ADJFACTOR=0.99;
    
    register unsigned int i,j, N=Metric.rno();
    unsigned int Viol=0;
    register double Sqroots, Dgshf=0.0;

    // get minimal metric diag (or 0.0 if all non-neg)
    Dgshf=0.0; Viol=0;
    for (i=0; i<N; i++)
	if (Metric[i][i]<Dgshf)
	{
	    Dgshf=Metric[i][i];
	    Viol++;	/* correction needed */
	}

    // shift centredists if necessary
    if (Viol)
    {
	Dgshf*=2.0;
	for (i=0; i<N; i++)
	    Metric[i][i]-=Dgshf;
    }
    Diagshf+=Dgshf;	// collect the total shift
    
    // generate, check and adjust off-diagonals
    for (i=1; i<N; i++)
	for (j=0; j<i; j++)	// Cosine rule
	{
	    Sqroots=sqrt(Metric[i][i]*Metric[j][j]);
	    if (Metric[i][j]<-Sqroots)
	    {
		Viol++;
		Metric[i][j]=-ADJFACTOR*Sqroots;
	    }
	    else if (Metric[i][j]>Sqroots)
	    {
		Viol++;
		Metric[i][j]=ADJFACTOR*Sqroots;
	    }
	}
    return(Viol);
}
// END of trieq_bal()

// ---- Metric matrix conversions ----

/* dist_metric(): calculates the metric matrix Metric from a
 * matrix of squared interpoint distances Dist. The 2-parameter
 * "standalone" version adjusts the size of Metric if necessary
 * and calculates the squared distances from the centroid internally.
 * The 3-parameter version (called by trineq_filter() only) does
 * not perform dist checks and expects a prefabricated Cdist2
 * squared point-centroid distance vector as input. Private
 */
void Iproj_::dist_metric(const Trimat_& Dist, Trimat_& Metric)
{
    register unsigned int i,j;
    unsigned int N=Dist.rno();
    Vector_ Cdist2(N);
    
    Metric.set_size(N);	// adjust if necessary
    centre_dist(Dist, Cdist2);	// get distances from centroid
    
    Metric[0][0]=Cdist2[0];	// same as 3-parameter version from here on
    for (i=1; i<N; i++)
    {
	Metric[i][i]=Cdist2[i];
	for (j=0; j<i; j++)
	    Metric[i][j]=(Cdist2[i]+Cdist2[j]-Dist[i][j])/2.0;
    }	
}

void Iproj_::dist_metric(const Trimat_& Dist, const Vector_& Cdist2, Trimat_& Metric)
{
    register unsigned int i,j;
    unsigned int N=Dist.rno();
    
    Metric[0][0]=Cdist2[0];
    for (i=1; i<N; i++)
    {
	Metric[i][i]=Cdist2[i];
	for (j=0; j<i; j++)
	    Metric[i][j]=(Cdist2[i]+Cdist2[j]-Dist[i][j])/2.0;
    }	
}
// END of dist_metric()

/* metric_dist(): the inverse of dist_metric() above.  Private */
void Iproj_::metric_dist(const Trimat_& Metric, Trimat_& Dist)
{
    register unsigned int i,j;
    unsigned int N=Metric.rno();
    
    Dist.set_size(N);
    for (i=0; i<N; i++)
    {
	Dist[i][i]=0.0;
	for (j=0; j<i; j++)
	    Dist[i][j]=Metric[i][i]+Metric[j][j]-2.0*Metric[i][j];
    }
}
// END of metric_dist()

/* centre_dist(): calculates the squared distances of points from their
 * common centroid, given the squared distance matrix Dist.
 * All points have the same mass.
 * Based on Lagrange's Theorem. No dim checks.
 * Return value: the no. of Cdist2 elements <0.0 (indicates non-metricity). Private
 */
unsigned int Iproj_::centre_dist(const Trimat_& Dist, Vector_& Cdist2)
{
    register double Trisum,Isum;
    register unsigned int i,j,k;
    unsigned int Cderr=0, N=Dist.rno();

    // get lower triangle squared sum
    Trisum=0.0;
    for (j=0; j<N; j++)
	for (k=0; k<j; k++)
	    Trisum+=Dist[j][k];
    Trisum/=(N*N);

    // get squared distance sums for the i-th point
    for (i=0; i<N; i++)
    {
	Isum=0.0;
	for (j=0; j<i; j++)
	    Isum+=Dist[i][j];
	for (j=i+1; j<N; j++)
	    Isum+=Dist[j][i];
	Cdist2[i]=Isum/N-Trisum;
	if (Cdist2[i]<0.0)
	{
	    cerr<<"? centre_dist(): Cdist2["<<i<<"]="<<Cdist2[i]<<endl;
	    Cderr++;
	}
    }
    return(Cderr);
}
// END of centre_dist()

/* centre_dist(): same as above but the masses of the points are supplied
 * in a vector Mass. Based on the "weighed" Lagrange's Theorem.
 */
unsigned int Iproj_::centre_dist(const Trimat_& Dist, const double *Mass, Vector_& Cdist2)
{
    register double Trisum, Isum, Msum;
    register unsigned int i,j,k;
    unsigned int Cderr=0, N=Dist.rno();

    // get lower triangle squared sum and sum of masses
    Trisum=Msum=0.0;
    for (j=0; j<N; j++)
    {
	Msum+=Mass[j];
	for (k=0; k<j; k++)
	    Trisum+=Mass[j]*Mass[k]*Dist[j][k];
    }
    Trisum/=(Msum*Msum);

    // get squared distance sums for the i-th point
    for (i=0; i<N; i++)
    {
	Isum=0.0;
	for (j=0; j<i; j++)
	    Isum+=Mass[j]*Dist[i][j];
	for (j=i+1; j<N; j++)
	    Isum+=Mass[j]*Dist[j][i];
	Cdist2[i]=Isum/Msum-Trisum;
	if (Cdist2[i]<0.0)
	{
	    cerr<<"? centre_dist(): Cdist2["<<i<<"]="<<Cdist2[i]<<endl;
	    Cderr++;
	}
    }
    return(Cderr);
}
// END of centre_dist()

// ---- Auxiliaries ----

/* sub_matrix(): selects a submatrix from the matrix Mat and copies
 * the elements to Submat. All entries between the active points
 * in Act are copied and the size of Submat is adjusted dynamically.
 * Return value: the no. of active points. Private
 */
unsigned int Iproj_::sub_matrix(const Trimat_& Mat, const Bits_& Act, 
	Trimat_& Submat)
{
    register unsigned int i, j, Rno=Act.on_no();
    register int di, dj;
    Submat.set_size(Rno);
    
    for (di=i=0; i<Mat.rno(); i++, di++)
    {
	if (!Act.get_bit(i)) { di--; continue; }	// inactive point, skip
	for (dj=j=0; j<=i; j++, dj++)
	{
	    if (!Act.get_bit(j)) { dj--; continue; }	// inactive
	    Submat[di][dj]=Mat[i][j];	// copy matrix entries
	}
    }
    return(Rno);
}
// END of sub_matrix()

/* make_locdist(): calculates the interpoint distances among the
 * active points in Locals and puts these into Locdist. 
 * Locdist is assumed to have been zeroed before the call. Private
 */
void Iproj_::make_locdist()
{
    register unsigned int i, j;
    register int di, dj;
    
    for (di=i=0; i<Rno; i++, di++)
    {
	if (!Locals.active(i)) { di--; continue; }	// inactive point, skip
	for (dj=j=0; j<=i; j++, dj++)
	{
	    if (!Locals.active(j)) { dj--; continue; }	// inactive
	    Locdist[i][j]=diff_len2(Locals[di], Locals[dj]);	// calc distances
	}
    }
}
// END of make_locdist() 

/* apply_locdist(): if Locdist[i][j]!=0.0, then the entry in
 * Dist[i][j] is replaced by it. Use for restoring the Euclidean
 * intra-cluster distances during triangle inequality smoothing
 * prior to skeleton embedding. Assumes that Dist is Rno x Rno. Private
 */
void Iproj_::apply_locdist(Trimat_& Dist) const
{
    register unsigned int i, j;
    register double Ld;
    
    for (i=0; i<Rno; i++)
	for (j=0; j<i; j++)
	{
	    Ld=Locdist[i][j];
	    if (Ld>0.0) Dist[i][j]=Ld;
	}
}
// END of apply_locdist() 

#undef DEBUG
#undef SKELSAVE

// ==== END OF METHODS Iproj.c++ ====
