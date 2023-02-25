#ifdef USE_OPENGL_GRAPHICS

/* ==== FUNCTIONS cadraw.c ==== */

/* ANSI C + X11 + OpenGL, 15-Oct-1997, Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <math.h>
#include <GL/glu.h>

/* ---- MODULE HEADER ---- */

#include "cadraw.h"

/* NOTE: SGI provides single-precision floating point functions
 * such as sqrtf() etc. Some machines (SUNs in particular) don't
 * know about this. Get around by the following macro
 */
#ifdef NO_MATHFLOATFUNC
#define sqrtf sqrt
#define fabsf fabs
#endif

/* ==== FUNCTIONS ==== */

/* ---- Data initialisation ---- */

/* create_drawchain(): creates a Drawchain_ object,
 * allocates storage for N points inside, 
 * and sets the coordinate minima and maxima to very high
 * and very low values.
 * Return value: the new object (NULL on failure).
 */
Drawchain_ *create_drawchain(int N)
{
    Drawchain_ *D=NULL;
    int i;
    
    N=fabs(N);
    if (!N) return(NULL);
    if (NULL==(D=(Drawchain_*) malloc(sizeof(Drawchain_))) ||
	NULL==(D->Coords=(Coord_ *) calloc(N, sizeof(Coord_))))
    {
	fputs("\n? init_drawchain(): Cannot allocate\n", stderr);
	return(NULL);
    }
    D->Cono=N;
    reset_drawlimits(D);
    return(D);
}
/* END of create_drawchain() */

/* delete_drawchain(): frees up the memory associated with *Drawchain
 * and removes the object itself.
 */
void delete_drawchain(Drawchain_ *Drawchain)
{
    if (Drawchain!=NULL)
    {
	free(Drawchain->Coords);
	free(Drawchain);
    }
}
/* END of delete_drawchain() */

/* reset_drawlimits(): clears the minimum and maximum coordinate
 * values so that a new minimum/maximum determination could be run.
 */
void reset_drawlimits(Drawchain_ *D)
{
    int i;
    if (D==NULL) return;
    for (i=0; i<3; i++)
    {
	D->Mincoord[i]=1e38;
	D->Maxcoord[i]=-1e38;
    }
}
/* END of reset_drawlimits() */

/* calc_drawlimits(): calculates the coordinate minima and maxima. */
void calc_drawlimits(Drawchain_ *D)
{
    register unsigned int i;
    register GLfloat X, Y, Z;
    
    if (D==NULL) return;
    
    reset_drawlimits(D);
    for (i=0; i<D->Cono; i++)
    {
	X=D->Coords[i].X[0];
	if (X<D->Mincoord[0]) D->Mincoord[0]=X;
	if (X>D->Maxcoord[0]) D->Maxcoord[0]=X;
	Y=D->Coords[i].X[1];
	if (Y<D->Mincoord[1]) D->Mincoord[1]=Y;
	if (Y>D->Maxcoord[1]) D->Maxcoord[1]=Y;
	Z=D->Coords[i].X[2];
	if (Z<D->Mincoord[2]) D->Mincoord[2]=Z;
	if (Z>D->Maxcoord[2]) D->Maxcoord[2]=Z;
    }
}
/* END of calc_drawlimits() */

/* ---- Graphics initialisation ---- */

/* init_cadraw(): sets up the GL context in *Winfo so that
 * z-buffering and antialiasing (with a line width Linewidth)
 * is enabled.
 */
void init_cadraw(Windowinfo_ *Winfo, GLfloat Linewidth)
{
    /* connect the context to the window */
    glXMakeCurrent(Winfo->Dpy, Winfo->Win, Winfo->Ctx);
    
    /* set up RGBA line antialiasing */
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glLineWidth(Linewidth);
    
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);	/* z-buffer */
    glClearColor(0.0, 0.0, 0.0, 0.0);	/* clear colour is pitch black */
}
/* END of init_cadraw() */

/* set_perspective(): calculate a bounding box for the chain
 * in *Drawchain, using the minimal and maximal coordinate
 * data therein. 
 * The viewing angle is specified in Viewangle in _degrees_ .
 * The matrix mode is set to GL_MODELVIEW upon return.
 */
void set_perspective(Windowinfo_ *Winfo, const Drawchain_ *Drawchain, 
	GLfloat Viewangle)
{
    static const float DEG_TO_RAD=0.017453293;
    
    static GLfloat Boxcenter[3];
    unsigned int Width, Height;
    XWindowAttributes Winattr;
    GLfloat Aspect;
    GLdouble Radius, Viewdist, Temp,
	Tangent=tan(0.5*DEG_TO_RAD*Viewangle);
    int i;
    
    /* calculate the bounding box center coords */
    for (i=0; i<3; i++)
	Boxcenter[i]=0.5*(Drawchain->Maxcoord[i]+Drawchain->Mincoord[i]);
    
    /* calculate the bounding sphere radius */
    Radius=0.0;
    for (i=0; i<3; i++)
    {
	Temp=Drawchain->Maxcoord[i]-Boxcenter[i];
	Temp*=Temp;
	Radius+=Temp;
    }
    Radius=sqrtf(Radius);
    Viewdist=Radius/Tangent;
    
    /* connect the context to the window */
    glXMakeCurrent(Winfo->Dpy, Winfo->Win, Winfo->Ctx);
    
    /* define the perspective */
    XGetWindowAttributes(Winfo->Dpy, Winfo->Win, &Winattr);
    Width=Winattr.width;
    Height=Winattr.height;
    Aspect=Width/Height;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-Radius, Radius, -Radius, Radius, Viewdist-Radius, Viewdist+Radius);
    
    /* shift the scene: X,Y to origin, Z to negative view distance */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-Boxcenter[0], -Boxcenter[1], -Boxcenter[2]-Viewdist);
    glViewport(0, 0, Width, Height);
}
/* END of set_perspective() */

/* ---- Drawing ---- */

/* draw_chain(): draws the C-alpha chain in *Drawchain
 * into the window *Winfo, using 3D perspective projection.
 * no rotations, antialiased line width of Linewidth.
 */
void draw_chain(Windowinfo_ *Winfo, const Drawchain_ *Drawchain)
{
    int i;
    
    /* connect the context to the window */
    glXMakeCurrent(Winfo->Dpy, Winfo->Win, Winfo->Ctx);
    
    /* drawing into back buffer if double-buffer mode is on */
    glDrawBuffer((Winfo->Dblbuffer)? GL_BACK: GL_FRONT);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    
    /* draw a polyline */
    glBegin(GL_LINE_STRIP);
    for (i=0; i<Drawchain->Cono; i++)
    {
	glColor3fv(Drawchain->Coords[i].Col);
	glVertex3fv(Drawchain->Coords[i].X);
    }
    glEnd();
    
    /* flush graphics */
    if (Winfo->Dblbuffer)
	glXSwapBuffers(Winfo->Dpy, Winfo->Win);
    else glFlush();
}
/* END of draw_chain() */

/* ==== END OF FUNCTIONS cadraw.c ==== */

#endif /* USE_OPENGL_GRAPHICS */
