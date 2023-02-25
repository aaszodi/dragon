#ifdef USE_OPENGL_GRAPHICS

#ifndef GRAPHICS_CLASS
#define GRAPHICS_CLASS

// ==== PROJECT DRAGON: HEADER Graphics.h ====

/* Graphics output routines which interface to the C/OpenGL routines
 * in "matplot" and "cadraw".
 */

// SGI C++ 7.1 + OpenGL, IRIX 6.2, 17-May-1998. Andris Aszodi

// ---- STANDARD HEADERS ----

#include <stdlib.h>

// ---- MODULE HEADERS ----

#include "Polymer.h"
#include "cadraw.h"
#include "matplot.h"

// ---- UTILITY HEADERS ----

#include "Points.h"

// ==== CLASSES ====

/* Class Graphics_: a very simple class that holds a ptr to a 
 * C struct called Drawchain_ (cf. "cadraw.h"). It can be asked to
 * draw a 3D C-alpha structure on the screen of an SGI machine.
 * It also plots distance matrices (cf. "matplot.h").
 */
class Graphics_
{
    // data
    private:
    
    static const double MINDIST, MAXDIST;	// colour distance limits
    static const int DXORIG, DYORIG, EXORIG, EYORIG, 
	    MXORIG, MYORIG, MATSIZE, MOVIESIZE;	// win positions and sizes
    
    Drawmatrix_ *Dmat, *Emat;  // "distance space" and "euclidean space" matrices
    Drawchain_ *Movie;	// the draw structure ptr
    unsigned int Size;	// residue no.
    Windowinfo_ Distwin, Euclwin, Moviewin;	// distmat and structure window IDs
    bool DisplayOK;	// false if one of the windows could not be opened
    
    // methods
    public:
    
	// constructors
    /* Init to empty (default). Copy and assignment are disabled. */
    Graphics_(): Dmat(NULL), Emat(NULL), Movie(NULL),
	    Size(0), DisplayOK(true)
    {
	Distwin.Dpy=Euclwin.Dpy=Moviewin.Dpy=NULL;  // no windows yet
    }
        
	// destructor
    ~Graphics_()
    {
	close_window();
	delete_drawmat(Dmat); Dmat=NULL;
	delete_drawmat(Emat); Emat=NULL;
	delete_drawchain(Movie); Movie=NULL;
    }
    
	// update and display
    /* update_polymer(): updates the calling object so that it can
     * display a structure specified by Polymer (as regards length 
     * and phobicity info). 
     */
    void update_polymer(const Polymer_& Polymer);
    
    /* display_dist(), display_eucl(): draw the Distance Space or
     * Euclidean Space distance matrices, respectively. The windows
     * are opened if necessary.
     */
    void display_dist(const Trimat_& Distmat);
    void display_eucl(const Trimat_& Distmat);

    /* display_coords(): copies the Xyz coordinates into Movie's
     * appropriate array for display, provided Xyz is 3-dimensional.
     * Opens a graphics window and draws the structure.
     * If not in 3D, then no action is taken.
     */
    void display_coords(const Points_& Xyz);
    
    /* close_window(): closes the windows if they were open. */
    void close_window();
    
    // Old-style matrix routines
    private:
    void copy_mat(const Trimat_& Newmat, Drawmatrix_ *Oldmat) const;

    // "forbidden methods"
    private:
    Graphics_(const Graphics_&);
    Graphics_& operator=(const Graphics_&);
};
// END OF CLASS Graphics_

// ==== END OF HEADER Graphics.h ====

#endif	/* GRAPHICS_CLASS */
#endif	/* USE_OPENGL_GRAPHICS */
