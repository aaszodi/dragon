#ifndef DSLCLU_HEADER
#define DSLCLU_HEADER

/* ==== HEADER dslclu.h ==== */

/* Single-linkage clustering of points with metric distances. */

/* ANSI C, IRIX 5.3, 14. Mar. 1996. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

/* ---- MODULE HEADERS ---- */

#include "matrix.h"

/* ---- TYPEDEFS ---- */

/* Dslclu_: the cluster object. Holds an array of int-s
 * which index the clustered things, pointers to two
 * sub-clusters (both are NULL if this is a leaf, single
 * sub-clusters not allowed) and the metric distance between
 * the sub-clusters.
 */
typedef struct dslclu_
{
    int *Members;   /* index of things in this cluster */
    int No;	    /* number of things */
    struct dslclu_ *Sub1, *Sub2;    /* ptrs to sub-clusters */
    float Dist;	    /* distance between sub-clusters */
} Dslclu_;

/* ---- PROTOTYPES ---- */

/* make_dslclus(): constructs a tree of Dslclu_ objects,
 * given the metric distance matrix Dist and Thingno, the number of
 * things (indexed between 0 and Thingno-1). The metricity
 * of Dist is assumed but NOT tested. Dist is assumed to
 * be Thingno x Thingno but there's no way of testing this :-)
 * The new distances between a freshly merged cluster and everybody else is
 * calculated as if the new cluster were a weighted average
 * of its subclusters.
 * Return value: a ptr to the root of the tree or NULL on error.
 */
Dslclu_ *make_dslclus(Trimat_ Dist, int Thingno);

/* clu_remove(): deletes the Members array of the cluster
 * pointed to by Clu. Recursively calls itself for the subclusters
 * as well if they exist and removes them, too. However, 
 * the object pointed to by Clu itself is NOT removed.
 */
void clu_remove(Dslclu_ *Clu);

/* print_dslclus(): prints the whole cluster tree with root pointed
 * to by Clu to Out.
 */
void print_dslclus(const Dslclu_ *Clu, FILE *Out);

/* ==== END OF HEADER dslclu.h ==== */

#endif	/* DSLCLU_HEADER */
