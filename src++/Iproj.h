#ifndef IPROJ_CLASS
#define IPROJ_CLASS

// ==== PROJECT DRAGON: HEADER Iproj.h ====

/* The Hierarchic Inertial Projection. */

// SGI C++ 4.0, IRIX 5.3, 19. Jul. 1996. Andris Aszodi 

// ---- STANDARD HEADERS ----

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>

// ---- UTILITY HEADERS ----

#include "Array.h"
#include "Vector.h"
#include "Matrix.h"
#include "Trimat.h"
#include "Points.h"

// ==== CLASSES ====

/* Class Iproj_ : the Inertial Projection class. Stores the local
 * coordinates and local moments of inertia for the clusters in
 * the point set to be projected. Can be asked to perform a
 * "full projection" when a complete distance matrix is
 * projected, or a "skeleton projection" afterwards when the 
 * local structures are translated and rotated as rigid bodies.
 */
class Iproj_
{
    // data
    private:
    
    Points_ Locals, Imoms;  // local cluster coordinates and moments of inertia
    Bits_ *Clusters;	    // the cluster membership
    unsigned int *Ptclu, *Ptoffs, *Cluoffs; // idx offset arrays, cf. make_offsets()
    unsigned int Rno, Cluno;	// residue and cluster numbers
    Trimat_ Skmet, Locdist;	// skeleton metric matrix and local Euclidean distances
    Points_ Skxyz;	// skeleton coordinates
    unsigned int Maxlocdim, Sksize; // maximal local dimension, skeleton size
    double Diagshf;	// diagonal shift
    
    // methods
    public:
    
	// constructors
    /* Inits to perform projections on a point set made up of Resno points */
    Iproj_(unsigned int Resno);
    
	// destructor
    ~Iproj_()
    {
	delete [] Clusters; delete [] Ptclu; 
	delete [] Ptoffs; delete [] Cluoffs;
    }
    
	// size
    /* set_size(): sets the calling object to Resno residues.
     * Works much in the same way as the corresponding constructor.
     * If Cl is 0 (the default), then the cluster number is calculated
     * within; any other number is interpreted as the required number
     * of clusters.
     */
    void set_size(unsigned int Resno, unsigned int Cl=0);
    
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
    unsigned int make_clusters(unsigned int Cno=0);
    unsigned int make_clusters(const Array_<Bits_>& Clus);

	// access
    /* ### For debug only. */
    const Bits_& clusters(unsigned int ci) const { return(Clusters[ci]); }
    unsigned int cluno() const { return(Cluno); }
    
	// projections
    
    /* full_project(): performs the Hierarchic Inertial Projection
     * on a point set. Dist holds the squared interpoint distances.
     * The projections will use an Evfract-th fraction of
     * the sum of all positive eigenvalues and will project into a
     * Dim<Oldim-dimensional Euclidean space. The Cartesian coordinates
     * will be put into Xyz (the sizes of which is set within).
     * Return value: the new dimension.
     */
    unsigned int full_project(Trimat_& Dist,
	    double Evfract, unsigned int Oldim, Points_& Xyz);
    
    /* cluster_project(): given a distance matrix Dist,
     * the local cluster coordinates and moments of
     * inertia are updated. Also sets the maximal local dimension
     * and the skeleton size.
     * Dist is modified: after local projection, the new local distances
     * (which are metric now) are put back into it.
     * No action is taken when invoked on 1-cluster sets.
     */
    void cluster_project(Trimat_& Dist, unsigned int Oldim);
    
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
    unsigned int skel_project(Trimat_& Dm, 
	    double Evfract, unsigned int Oldim, Points_& Xyz);
    
    // private methods
    private:
	// size
    void make_offsets();
    
	// reconstruction
    void flesh_skel(const Trimat_& Dist, Points_& Xyz);
    static float clu_qual(const Bits_& Clu, const Points_& Xyz, const Trimat_& Dist);

	// projections
    unsigned int metric_project(const Trimat_& Metric,  
	    double Evfract, unsigned int Mindim, unsigned int Oldim, 
	    Points_& Xyz, Vector_ *Moms=NULL);
    
	// scalar products
    void make_skmet(const Trimat_& Metric);
    void ctr_prod(const Trimat_& Metric, Matrix_& Abprods) const;
    double iv_ctrprod(const Matrix_& Abprods,
	    unsigned int Ctridx, unsigned int Coffs, unsigned int Coord) const;
    void aibj_prod(const Trimat_& Metric, const Matrix_& Abprods, 
	    unsigned int Aidx, unsigned int Bidx, 
	    unsigned int Aoffs, unsigned int Boffs, 
	    Matrix_& Aibj) const;

	// triangle inequality balance
    unsigned int trineq_filter(Trimat_& Dist, Trimat_& Metric,
	unsigned int Tsmcyc=0, const double *Mass=NULL);
    unsigned int trieq_bal(Trimat_& Metric);

	// metric matrix conversions
    static void dist_metric(const Trimat_& Dist, Trimat_& Metric);
    static void dist_metric(const Trimat_& Dist, const Vector_& Cdist2, Trimat_& Metric);
    static void metric_dist(const Trimat_& Metric, Trimat_& Dist);
    static unsigned int centre_dist(const Trimat_& Dist, Vector_& Cdist2);
    static unsigned int centre_dist(const Trimat_& Dist, const double *Mass, Vector_& Cdist2);
	
	// auxiliaries
    static unsigned int sub_matrix(const Trimat_& Mat, const Bits_& Act, 
	Trimat_& Submat);
    
    void make_locdist();
    void apply_locdist(Trimat_& Dist) const;

    // "forbidden methods"
    private:
    
    Iproj_();
    Iproj_(const Iproj_&);
    Iproj_& operator=(const Iproj_&);
    
};
// END OF CLASS Iproj_;

// ==== END OF HEADER Iproj.h ====
#endif	/* IPROJ_CLASS */
