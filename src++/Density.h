#ifndef __DENSITY_HEADER__
#define __DENSITY_HEADER__

// ==== PROJECT DRAGON: HEADER Density.h ====

/* Adjusts molecular density. */

// SGI C++ 4.0, IRIX 5.3, 4. July 1995. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- MODULE HEADERS ----

#include "Points.h"
#include "Pieces.h"

// ---- UTILITY HEADERS ----

#include "Trimat.h"

// ---- PROTOTYPES ----

/* scale_distdens(): uses Willie's "moment scaling" to adjust the
 * density of the distance matrix to the density expected for a
 * 3D spherical protein. He has calculated the first and second moments
 * of the distance distribution between any pairs of points within
 * a solid sphere of radius Rmax. Here the moments are found and a
 * scale factor is calculated as the average of the ratios of the
 * expected moments to the observed ones.
 * Return value: the scaling factor.
 */
double scale_distdens(Trimat_& Dist, double Rmax);

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
double proj_dens(const Trimat_& Dist, const Pieces_& Pieces, Points_& Xyz);

/* ellips_dens(): adjusts the density of the 3D Euclidean point set Xyz
 * to match the expected density Expdens. Fits an ellipsoid to the points
 * so that it contains 90% of them,  works out an adjustment factor and
 * updates the points. Single-cluster point sets are isotropically adjusted, 
 * multicluster sets are updated by moving the clusters as rigid bodies.
 * Return value: the adjustment factor. If the dimension is not 3D, 
 * then no action is taken and 0.0 is returned.
 */
double ellips_dens(double Expdens, const Pieces_& Pieces, Points_& Xyz);

// ==== END OF HEADER Density.h ====
#endif
