#ifdef USE_OPENGL_GRAPHICS

#ifndef CADRAW_H
#define CADRAW_H

/* ==== PROJECT DRAGON: HEADER cadraw.h ==== */

/* Draws a colour-coded C-alpha polypeptide chain.
 * Only static images are supported -- no rotation.
 */

/* ANSI C + X11 + OpenGL, 8-Oct-1997, Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glx.h>

/* ---- MODULE HEADERS ---- */

#include "glxwinutils.h"

/* ---- TYPEDEFS ---- */

/* Coord_: stores the C-alpha atom coordinates
 * and the corresponding colour to be drawn.
 */
typedef struct
{
    GLfloat X[3];	/* coordinates */
    GLfloat Col[3];	/* RGB colour */
}
Coord_ ;

/* Drawchain_: the whole chain to be drawn.
 * Holds an array of Coord_ items and its length.
 * When the coordinates are initialised, calculate
 * the minimal and maximal values and store them in Mincoord
 * and Maxcoord.
 */
typedef struct
{
    Coord_ *Coords;	/* atom coordinates & colour */
    int Cono;		/* number of points */
    GLfloat Mincoord[3], Maxcoord[3];	/* min and max values for X,Y,Z */
}
Drawchain_ ;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Data initialisation ---- */

/* create_drawchain(): creates a Drawchain_ object,
 * allocates storage for N points inside, 
 * and sets the coordinate minima and maxima to very high
 * and very low values.
 * Return value: the new object (NULL on failure).
 */
Drawchain_ *create_drawchain(int N);

/* delete_drawchain(): frees up the memory associated with *Drawchain
 * and removes the object itself.
 */
void delete_drawchain(Drawchain_ *Drawchain);

/* reset_drawlimits(): clears the minimum and maximum coordinate
 * values so that a new minimum/maximum determination could be run.
 */
void reset_drawlimits(Drawchain_ *D);

/* calc_drawlimits(): calculates the coordinate minima and maxima. */
void calc_drawlimits(Drawchain_ *D);

/* ---- Graphics initialisation ---- */

/* init_cadraw(): sets up the GL context in *Winfo so that
 * z-buffering and antialiasing (with a line width Linewidth)
 * is enabled.
 */
void init_cadraw(Windowinfo_ *Winfo, GLfloat Linewidth);

/* set_perspective(): calculate a bounding box for the chain
 * in *Drawchain, using the minimal and maximal coordinate
 * data therein. 
 * The viewing angle is specified in Viewangle in _radians_ .
 * The matrix mode is set to GL_MODELVIEW upon return.
 */
void set_perspective(Windowinfo_ *Winfo, const Drawchain_ *Drawchain, 
	GLfloat Viewangle);

/* ---- Drawing ---- */

/* draw_chain(): draws the C-alpha chain in *Drawchain
 * into the window *Winfo, using 3D perspective projection.
 * no rotations, antialiased line width of Linewidth.
 */
void draw_chain(Windowinfo_ *Winfo, const Drawchain_ *Drawchain);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER cadraw.h ==== */
#endif	/* CADRAW_H */

#endif	/* USE_OPENGL_GRAPHICS */
