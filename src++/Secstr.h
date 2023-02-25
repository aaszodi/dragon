#ifndef __SECSTR_CLASSES__
#define __SECSTR_CLASSES__

// ==== PROJECT DRAGON: HEADER Secstr.h ====

/* Classes for handling secondary structures: H-bond topology
 * and ideal geometry. Chain topology comes from the Segment module.
 * This header should be included into all modules using secondary
 * structure features: it includes, in turn, all necessary secstr
 * library headers.
 */

/* The model chains are divided into segments which correspond to the
 * secondary structure layout in the current implementation. The classes
 * in the "Segment" module represent the topology of the secstr segments, while
 * the classes in this module represent the geometry (ideal
 * distances and structure w/ chirality). 
 *
 * The inheritance graph of the whole family is as follows:-
 *
 *	         [Segmbase_]
 *                    |
 *	    +---------+----------+
 *	    :         :          :
 *	    V         V          V
 *	Linsegm_  [Sstrbase_]  Sheet_
 *         |           |          |
 *         +--------+  +--------+ |
 *         |        |  |        | |
 *	   V        V  V        V V
 *	Strand_     Helix_      Beta_
 * 
 * ( [] enclose abstract base classes, : denotes virtual ancestry)
 */

// SGI C++ 4.0, IRIX 5.3, 11. Aug. 1995. (C) Andras Aszodi

// ---- LIBRARY HEADERS ----

#include "Segment.h"
#include "Sstrbase.h"
#include "Helix.h"
#include "Beta.h"

// ==== END OF HEADER Secstr.h ====
#endif
