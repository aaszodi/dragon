#ifdef USE_OPENGL_GRAPHICS

/* ==== PROJECT DRAGON: HEADER matplot.h ==== */

#ifndef MATPLOT_H
#define MATPLOT_H

/* ==== HEADER matplot.h ==== */

/* OpenGL routines for visualising general rectangular MxN matrices. */

/* ANSI C+ OpenGL, 18-May-1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glx.h>

/* ---- MODULE HEADERS ---- */

#include "glxwinutils.h"

/* ---- TYPEDEFS ---- */

/* Drawmatrix_: this structure stores a ptr to the matrix
 * to be drawn, its size and the lower and upper limits of
 * the values for which the colour coding must be applied.
 */
typedef struct
{
    double **Mat;   /* the data structure */
    unsigned Row, Col;	/* row and column size */
    double Lowval, Upval; /* lowest and highest matrix element */
    GLboolean Resizeneeded;     /* GL_TRUE if the display window was resized */
}
Drawmatrix_;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* create_drawmat(): creates a Drawmatrix_ structure
 * which has a size RxC and the lower and upper limits
 * are set to Low and Up, respectively.
 * Returns a ptr to the newly allocated object.
 */
Drawmatrix_ *create_drawmat(unsigned int R, unsigned int C, 
	double Low, double Up);

/* delete_drawmat(): frees up the storage associated with
 * *Drawmat, then removes the object itself.
 */
void delete_drawmat(Drawmatrix_* Drawmat);

/* plot_mat: creates a simple colour-coded dot representation of the
 * matrix stored in the Drawmat structure.
 * Plotting is performed in the window identified by *Windowinfo.
 */
void plot_mat(Windowinfo_ *Winfo, Drawmatrix_ *Drawmat);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER matplot.h ==== */
#endif	/* MATPLOT_H */

#endif /* USE_OPENGL_GRAPHICS */
