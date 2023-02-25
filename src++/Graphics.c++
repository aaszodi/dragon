#ifdef USE_OPENGL_GRAPHICS

// ==== PROJECT DRAGON: METHODS Graphics.c++ ====

/* Graphics output routines which interface to the C/OpenGL routines
 * in "cadraw" and "matplot".
 */

// SGI C++ 7.1 + OpenGL, IRIX 6.2, 21-Oct-1997. Andris Aszodi

// ---- MODULE HEADER ----

#include "Graphics.h"
#include <math.h>
#include <string.h>

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this. Get around by the following macro
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#define fabsf fabs
#endif

// ---- STATIC INITIALISATION ----

const double Graphics_::MINDIST=3.5;
const double Graphics_::MAXDIST=15.0;

const int Graphics_::DXORIG=100;
const int Graphics_::DYORIG=100;
const int Graphics_::EXORIG=320;
const int Graphics_::EYORIG=100;
const int Graphics_::MXORIG=100;
const int Graphics_::MYORIG=340;
const int Graphics_::MATSIZE=200;
const int Graphics_::MOVIESIZE=400;

// ==== Graphics_ METHODS ====

// ---- Update and display ----

/* update_polymer(): updates the calling object so that it can
 * display a structure specified by Polymer (as regards length 
 * and actual phobicity info).
 */
void Graphics_::update_polymer(const Polymer_& Polymer)
{
    register unsigned int i, Ptno=Polymer.len()+2;   // room for N/C termini
    
    // change sizes
    if (Size!=Ptno)
    {
	// destroy old stuff
	delete_drawmat(Dmat); Dmat=NULL;
	delete_drawmat(Emat); Emat=NULL;
	delete_drawchain(Movie); Movie=NULL;
	Size=Ptno;

	// rebuild the matrices using the new Size
	Dmat=create_drawmat(Size, Size, MAXDIST, MINDIST);
	Emat=create_drawmat(Size, Size, MAXDIST, MINDIST);
	
	// construct the movie drawchain entry
	Movie=create_drawchain(Size);
    }
    
    // store colour information
    float Phob, Minphob=1e10, Maxphob=-1e10;
    for (i=0; i<Ptno-2; i++)	// get phobicity limits
    {
	Phob=Polymer.phob(i);
	if (Phob<Minphob) Minphob=Phob;
	if (Phob>Maxphob) Maxphob=Phob;
    }
    
    memcpy(Movie->Coords[0].Col, rainbow_ramp(Polymer.phob(0), Minphob, Maxphob), 3*sizeof(GLfloat));
    memcpy(Movie->Coords[Ptno-1].Col, rainbow_ramp(Polymer.phob(Ptno-3), Minphob, Maxphob), 3*sizeof(GLfloat));
    
    // store ordinary C-alpha atom colour codes
    for (i=1; i<Ptno-1; i++)
	memcpy(Movie->Coords[i].Col, rainbow_ramp(Polymer.phob(i-1), Minphob, Maxphob), 3*sizeof(GLfloat));
}
// END of update_polymer()

/* display_dist(), display_eucl(): draw the Distance Space or
 * Euclidean Space distance matrices, respectively. The windows
 * are opened if necessary.
 */
void Graphics_::display_dist(const Trimat_& Distmat)
{
    if (!DisplayOK) return;	// do not even try
    
    // shall the window be opened?
    if (Distwin.Dpy==NULL)
    {
	DisplayOK=create_glxwindow(&Distwin, DXORIG, DYORIG, MATSIZE, MATSIZE, "DRAGON:Dist");
	if (!DisplayOK)
	{
	    cerr<<"\n? Graphics_::display_dist(): No graphics\n";
	    return;
	}
    }
    copy_mat(Distmat, Dmat);	// copies square roots!
    plot_mat(&Distwin, Dmat);
}
// END of display_dist()

void Graphics_::display_eucl(const Trimat_& Distmat)
{
    if (!DisplayOK) return;	// do not even try
    
    // shall the window be opened?
    if (Euclwin.Dpy==NULL)
    {
	DisplayOK=create_glxwindow(&Euclwin, EXORIG, EYORIG, MATSIZE, MATSIZE, "DRAGON:Eucl");
	if (!DisplayOK)
	{
	    cerr<<"\n? Graphics_::display_eucl(): No graphics\n";
	    return;
	}
    }
    copy_mat(Distmat, Emat);	// copies square roots!
    plot_mat(&Euclwin, Emat);
}
// END of display_eucl()

/* display_coords(): copies the Xyz coordinates into Movie's
 * appropriate array for display, provided Xyz is 3-dimensional.
 * Opens a graphics window and draws the structure.
 * If not in 3D, then no action is taken.
 */
void Graphics_::display_coords(const Points_& Xyz)
{
    // are we in 3D, was graphics OK so far?
    if (Xyz.dim()!=3 || !DisplayOK) return;	// no action
    
    register unsigned int i;
    Coord_ *Coords=Movie->Coords;
    float *Cx;
    static GLenum Sizechanged=GL_FALSE;
    
    // copy coordinates
    for (i=0; i<Movie->Cono; i++)
    {
	Cx=Coords[i].X;
	Cx[0]=(float)Xyz[i][0];
	Cx[1]=(float)Xyz[i][1];
	Cx[2]=(float)Xyz[i][2];
    }
    
    // open the window
    if (Moviewin.Dpy==NULL)
    {
	DisplayOK=create_glxwindow(&Moviewin, MXORIG, MYORIG, MOVIESIZE, MOVIESIZE, "DRAGON:Movie");
	if (!DisplayOK)
	{
	    cerr<<"\n? Graphics_::display_coords(): No graphics\n";
	    return;
	}
	init_cadraw(&Moviewin, 3);  // linewidth etc.
	calc_drawlimits(Movie);    // prepare for perspective
	Sizechanged=GL_TRUE;	    // make a note to set perspective
    }
    
    // has the size changed?
    if (Sizechanged || 
	(read_events(&Moviewin) && Moviewin.Event.type==ConfigureNotify))
    {
	set_perspective(&Moviewin, Movie, 20.0);	// 20 deg window
	Sizechanged=GL_FALSE;
    }
    
    draw_chain(&Moviewin, Movie);	// draw the structure
}
// END of display_coords()

/* close_window(): closes the windows if they were open. */
void Graphics_::close_window()
{
    if (!DisplayOK) return;
    destroy_glxwindow(&Distwin);
    destroy_glxwindow(&Euclwin);
    destroy_glxwindow(&Moviewin);
}

// ---- Old matrix handling ----

void Graphics_::copy_mat(const Trimat_& Newmat, Drawmatrix_ *Oldmat) const
{
    if (Newmat.rno()!=Oldmat->Row || Newmat.rno()!=Oldmat->Col)
    {
	cerr<<"\n? Graphics_::copy_mat(): Dim mismatch: Oldmat=["
		<<Oldmat->Row<<","<<Oldmat->Col<<"], Newmat=["
		<<Newmat.rno()<<","<<Newmat.rno()<<"]\n";
	return;
    }
    
    register unsigned int i, j;
    for (i=0; i<Size; i++)
    {
	Oldmat->Mat[i][i]=sqrtf(fabsf(Newmat[i][i]));
	for (j=0; j<i; j++)
	    Oldmat->Mat[i][j]=Oldmat->Mat[j][i]=sqrtf(fabsf(Newmat[i][j]));
    }
}
// END of copy_mat()

// ==== END OF METHODS Graphics.c++ ====

#endif	/* USE_GRAPHICS */
