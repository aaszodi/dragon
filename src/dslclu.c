/* ==== dslclu ==== */

/* Single-linkage clustering of points with metric distances. */

/* ANSI C, IRIX 5.3, 14. Mar. 1996. Andris */

/* ---- HEADER ---- */

#include "dslclu.h"

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this: Define NO_MATHFLOATFUNC on the command line
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#endif

/* ---- PROTOTYPES ---- */

static float find_closest(Trimat_ Dist, int Size, int *Ci, int *Cj);
static void merge_clus(Dslclu_ *Clu1, Dslclu_ *Clu2, float Cludist);
static void update_distmat(Trimat_ Dist, int Size, int Ci, int Cj,
	float Wi, float Wj);
static int new_dists(Trimat_ Dist, int Size, int Ci, int Cj, 
	float Wi, float Wj, float Newdists[]);

static void init_dslclu(Dslclu_ *Clu, int Thing, int Thingno);
static Dslclu_ *clu_dup(const Dslclu_ *Clu);

/* ==== FUNCTIONS ==== */

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
Dslclu_ *make_dslclus(Trimat_ Dist, int Thingno)
{
    Dslclu_ *Clus=NULL, *Tree=NULL;
    Trimat_ Tempmat=NULL;
    int Cluno, Ci, Cj, i, j;
    float Closedist;
    
    if (Dist==NULL || Thingno<=0)
    {
	fprintf(stderr, "\n? make_dslclus(%p, %d): Silly parameter\n", 
	    (void*)Dist, Thingno);
	return(NULL);
    }
    
    /* create an initial cluster for each thing */
    Clus=(Dslclu_*) calloc(Thingno, sizeof(Dslclu_));
    for (i=0; i<Thingno; i++) 
	init_dslclu(Clus+i, i, Thingno);
    
    /* copy the original distmat to a temporary one which gets overwritten */
    Tempmat=alloc_trimat(Thingno);
    for (i=0; i<Thingno; i++)
	for (j=0; j<=i; j++)
	    Tempmat[i][j]=Dist[i][j];
    
    /* perform the actual clustering. In the end,
     * Clus[0] will be the root of the tree
     */
    for (Cluno=Thingno; Cluno>=2; Cluno--)
    {
	/* get pair to be merged */
	Closedist=find_closest(Dist, Cluno, &Ci, &Cj);

	/* update the array and the matrix,
	 * merge the pair: Clus[Cj] will be absorbed into Clus[Ci] 
	 */
	update_distmat(Dist, Cluno, Ci, Cj, Clus[Ci].No, Clus[Cj].No);
	merge_clus(Clus+Ci, Clus+Cj, Closedist);
	if (Cj<Cluno-1)    /* compact Clus[]: overlap/copy */
	    memmove(Clus+Cj, Clus+Cj+1, (Cluno-Cj-1)*sizeof(Dslclu_));
    }
    
    /* we're done: clean up */
    Tree=(Dslclu_*) malloc(sizeof(Dslclu_));
    *Tree=Clus[0];  /* root copied to Tree */
    free(Clus);
    free_matrix(Tempmat);
    return(Tree);
}
/* END of make_dslclus() */

/* ---- Clustering auxiliaries ---- */

/* find_closest(): return the indices in Ci and Cj of the smallest entry
 * in the Size x Size matrix Dist. I won't tell you which
 * entry is returned if there are more than one equally
 * smallest ones. The actual value is returned by the fn.
 */
static float find_closest(Trimat_ Dist, int Size, int *Ci, int *Cj)
{
    register int i, j;
    register float Smallest=FLT_MAX;
    
    for (i=0; i<Size; i++)
	for (j=0; j<i; j++)
	{
	    if (Dist[i][j]<Smallest)
	    {
		Smallest=Dist[i][j];
		*Ci=j; *Cj=i;	/* *Ci<=*Cj all the time */
	    }
	}
    return(Smallest);
}
/* END of find_closest() */

/* merge_clus(): the cluster pointed to by Clu2 will be merged into
 * Clu1. First both are duplicated, then Clu1 is updated so that
 * its two subclusters will be the duplicated copies and its
 * membership array will contain the union of those of the subclusters.
 * Their distance is Cludist which will be stored in the new Clu1.
 */
static void merge_clus(Dslclu_ *Clu1, Dslclu_ *Clu2, float Cludist)
{
    /* merge: here make sure non-terminals have 2 sub-nodes */
    Dslclu_ *Oldclu1=clu_dup(Clu1), *Oldclu2=clu_dup(Clu2);
    
    memcpy(Clu1->Members+Clu1->No, Clu2->Members, Clu2->No*sizeof(int));
    Clu1->No+=Clu2->No;
    Clu1->Dist=Cludist;
    Clu1->Sub1=Oldclu1; Clu1->Sub2=Oldclu2;	/* store subclusters */
    free(Clu2->Members); Clu2->Members=NULL;
}
/* END of merge_clus() */

/* update_distmat(): rewrites the Size x Size distance matrix Dist
 * after the Cj>Ci-th entry was merged into the Ci-th.
 * The Cj-th row (col) is eliminated altogether, and the distances
 * between the (new) Ci-th entry and the rest are updated.
 * This update is carried out so that the dists between the new
 * cluster and everybody else is calculated as if the new cluster
 * were the centroid of the two merged clusters. (Can be done
 * since the dist matrix is metric.) Nothing done if Ci>=Cj.
 * Wi and Wj are the combination weights of the two old clusters.
 */
static void update_distmat(Trimat_ Dist, int Size, int Ci, int Cj,
	float Wi, float Wj)
{
    register int i, j;
    int Isnew=0;
    float *Newdists=NULL;
    
    /* incomplete error checks */
    if (Ci>=Cj || Ci>=Size || Cj>=Size)
    {
	fprintf(stderr, "\n? update_distmat(..., %d, %d): Silly indices\n", 
		Ci, Cj);
	return;
    }
    
    /* calculate the distances from the new cluster to the rest:
     * if new_dists() returns 0 that means that Ci and Cj coincide
     * and no new dists had to be calculated
     */
    Newdists=(float *) calloc(Size-1, sizeof(float));
    Isnew=new_dists(Dist, Size, Ci, Cj, Wi, Wj, Newdists);
    
    /* move the section "below" the Cj-th row "up" by one */
    for (i=Cj; i<Size-1; i++)
	for (j=0; j<=i; j++)
	    Dist[i][j]=Dist[i+1][j];
    
    /* move the triangle "right" to the Cj-th col "to the left" by one */
    for (j=Cj; j<Size-1; j++)
    {
	for (i=j+1; i<Size; i++)
	    Dist[i][j]=Dist[i][j+1];
	Dist[j][j]=0.0;
    }
    
    /* copy the new distances into the Ci-th row/col */
    if (Isnew)	/* Ci, Cj didn't coincide */
    {
	for (j=0; j<=Ci; j++)
	    Dist[Ci][j]=Newdists[j];
	for (i=Ci+1; i<Size-1; i++)
	    Dist[i][Ci]=Newdists[i];
    }
    
    free(Newdists);
}
/* END of update_distmat() */

/* new_dists(): calculates the distances from the new Ci-th cluster
 * to everybody else and puts them into the Newdists[] array.
 * The distances are calculated so that the new cluster is a
 * point on the Ci:Cj segment, dividing it at a ratio of Wi:Wj.
 * The formula was obtained with the Cosine Rule (man's best friend).
 * On return, Newdists[Ci]==0.0, the Cj-th position is skipped
 * (in fact, Newdists is a Size-1 long array while Dist is still
 * the original Size x Size trimat).
 * Return value: 0 if the Ci:th and Cj:th points coincide
 * (in this case Newdists[] is not touched), 1 otherwise.
 */
static int new_dists(Trimat_ Dist, int Size, int Ci, int Cj, 
	float Wi, float Wj, float Newdists[])
{
    register float aa, bb, cc, dd, p1;
    register int i;
    
    /* check dist(ci,cj)==0 special case */
    aa=Dist[Cj][Ci];
    if (aa<1e-10)    /* assuming a>=0 */
	return(0);
    
    /* divide up "a" at Wi:Wj ratio */
    p1=Wj/(Wi+Wj); aa*=aa;
    aa*=(p1*p1-p1);
    
    /* calc new distances, assumes Cj>Ci */
    for (i=0; i<Ci; i++)   /* do the row */
    {
	bb=Dist[Ci][i]; bb*=bb;
	cc=Dist[Cj][i]; cc*=cc;
	dd=aa+bb*(1.0-p1)+cc*p1;
	Newdists[i]=sqrtf(dd);
    }
    Newdists[Ci]=0.0;	/* by def */
    for (i=Ci+1; i<Size; i++)	/* do the col */
    {
	if (i==Cj) continue;	/* skip Cj:th itself */
	bb=Dist[i][Ci]; bb*=bb;
	cc=(i<Cj)? Dist[Cj][i]: Dist[i][Cj]; cc*=cc;
	dd=aa+bb*(1.0-p1)+cc*p1;
	Newdists[(i<Cj)? i: i-1]=sqrtf(dd);
    }
    return(1);
}
/* END of new_dists() */

/* ---- Dslclu_ handling routines ---- */

/* init_dslclu(): sets up *Clu to store max. Thingno things,
 * and stores Thing as its only member. Would be a ctor in C++.
 */
static void init_dslclu(Dslclu_ *Clu, int Thing, int Thingno)
{
    Clu->Members=(int *) calloc(Thingno, sizeof(int));
    Clu->Members[0]=Thing;
    Clu->No=1; Clu->Dist=0.0;
    Clu->Sub1=Clu->Sub2=NULL;
}
/* END of init_dslclu() */

/* clu_dup(): duplicates the cluster pointed to by Clu and returns
 * a ptr to the new cluster. The Clu->Members array will be
 * duplicated and copied, the subclusters won't (only
 * the Sub1, 2 pointers).
 */
static Dslclu_ *clu_dup(const Dslclu_ *Clu)
{
    Dslclu_ *Newclu=(Dslclu_ *) malloc(sizeof(Dslclu_));
    *Newclu=*Clu;   /* use built-in ANSI memberwise copy */
    
    /* duplicate the member array */
    Newclu->Members=(int *) calloc(Newclu->No, sizeof(int));
    memcpy(Newclu->Members, Clu->Members, Clu->No*sizeof(int));
    return(Newclu);
}
/* END of clu_dup() */

/* clu_remove(): deletes the Members array of the cluster
 * pointed to by Clu. Recursively calls itself for the subclusters
 * as well if they exist and removes them, too. However, 
 * the object pointed to by Clu itself is NOT removed.
 */
void clu_remove(Dslclu_ *Clu)
{
    free(Clu->Members);	/* removes members */
    if (Clu->Sub1!=NULL)
    {
	clu_remove(Clu->Sub1);
	free(Clu->Sub1);    /* remove subcluster itself! */
    }
    if (Clu->Sub2!=NULL)
    {
	clu_remove(Clu->Sub2);
	free(Clu->Sub2);
    }
}
/* END -f clu_remove() */

/* print_dslclus(): prints the whole cluster tree with root pointed
 * to by Clu to Out.
 */
void print_dslclus(const Dslclu_ *Clu, FILE *Out)
{
    static const int NODEWIDTH=10;   /* "-[x.xe+00]" */
    static char *Hor=NULL, *Vert=NULL;    /* connecting lines etc. */
    static int Len=0;
    
    if (Clu==NULL) return;  /* bail out */
    
    /* init */
    if (!Len)
    {
	Hor=(char *) calloc(1, sizeof(char));	/* automatic \0-termination */
	Vert=(char *) calloc(1, sizeof(char));
    }
    
    /* leaf */
    if (Clu->Sub1==NULL && Clu->Sub2==NULL)
    {
	fprintf(Out, "%s-%d\n%s\n", Hor, Clu->Members[0], Vert);
	return;
    }
    
    /* non-leaf: always has two sub-nodes */
    if (Clu->Sub1!=NULL)    /* go to "left" */
    {
	Len+=NODEWIDTH;
	Hor=(char *) realloc(Hor, (Len+1)*sizeof(char));
	Vert=(char *) realloc(Vert, (Len+1)*sizeof(char));
	Hor[Len]=Vert[Len]='\0';
	sprintf(Hor+Len-NODEWIDTH, "-[%7.1e]", Clu->Dist);
	strcpy(Vert+Len-NODEWIDTH, "     |    ");
	print_dslclus(Clu->Sub1, Out);
    
	/* go to "right" */
	strcpy(Hor, Vert);
	strcpy(Hor+Len-NODEWIDTH, "     ^----");
	strcpy(Vert+Len-NODEWIDTH, "          ");
	print_dslclus(Clu->Sub2, Out);

	Len-=NODEWIDTH;	    /* shrink back to previous */
	if (Len)
	{
	    Hor=(char *) realloc(Hor, (Len+1)*sizeof(char));
	    Vert=(char *) realloc(Vert, (Len+1)*sizeof(char));
	    Hor[Len]=Vert[Len]='\0';
	}
	else	/* at topmost leaf */
	{
	    free(Hor); free(Vert);
	    Hor=Vert=NULL; Len=0;
	}
    }
    return;
}
/* END of print_dslclus() */

/* ==== END OF FUNCTIONS dslclu.c ==== */
