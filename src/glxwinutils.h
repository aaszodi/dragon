#ifdef USE_OPENGL_GRAPHICS

#ifndef GLXWINUTILS_H
#define GLXWINUTILS_H

/* ==== HEADER glxwinutils.h ==== */

/* A few utility routines to open a GL window under X. */

/* ANSI C + X11 + OpenGL, 17-Oct-1997. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glx.h>

/* ---- TYPEDEFS ---- */

/* Windowinfo_: this struct keeps all the magic
 * X/GL variables together which are needed to
 * describe a window into which we want to draw.
 */
typedef struct
{
    Display *Dpy;   /* X display */
    XVisualInfo *Visinfo;   /* X visual info */
    XSetWindowAttributes Winattr;   /* window attributes */
    Window Win;	    /* the window itself */
    GLXContext Ctx; /* drawing context */
    XEvent Event;   /* X events (resize, redraw...) */
    GLboolean Dblbuffer;    /* GL_TRUE if double buffering is OK */
} Windowinfo_;

/* Glxwucols_: these enums are symbolic constants for
 * frequently used RGB colour triplets. See glxwu_colour().
 * The last enum must always be GLXWU_COLNO.
 */
typedef enum {GLXWU_BLACK=0, GLXWU_BLUE, GLXWU_CYAN, 
	GLXWU_GREEN, GLXWU_YELLOW, GLXWU_RED, GLXWU_WHITE, 
	GLXWU_COLNO} Glxwucols_;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* create_glxwindow(): opens a GL/X window with its lower left
 * corner at Xorig:Yorig, size Xsize:Ysize, sets the title to Title and puts
 * all the gory details into *Windowinfo.
 * Aspect ratio is kept at Xsize/Ysize.
 * Returns GL_FALSE on error, GL_TRUE on success.
 */
GLenum create_glxwindow(Windowinfo_ *Winfo, int Xorig, int Yorig, 
	int Xsize, int Ysize, const char *Title);

/* destroy_glxwindow(): removes the window *Windowinfo from the screen. */
void destroy_glxwindow(Windowinfo_ *Winfo);

/* read_events(): reads the event queue associated with Winfo->Dpy "dry"
 * and returns GL_TRUE if Expose, VisibilityNotify or ConfigureNotify
 * events were detected. If only Expose was seen then Winfo->Event.type will
 * contain Expose, if at least one ConfigureNotify was seen then
 * Winfo->Event will contain the full description of it (to help
 * with resizing).
 */
GLenum read_events(Windowinfo_ *Winfo);

/* rainbow_ramp: converts a double value X to a colour in RGB
 * representation according to the following approx. scheme:
 * 
 * Colour  <--Black  Blue  Cyan   Green  Yellow  Red  White-->
 *                    |	    |       |       |      |
 * X (relative)       0 ...1/4.....1/2.....3/4.....1
 * X (absolute)      Lowval.....................Upval
 * Return value: const ptr to a 3-element GLfloat vector Colorvector
 * which is stored inside and is updated within calls.
 * Set the color with glColor3fv().
 * If Upval==Lowval, then black is returned.
 */
const GLfloat *rainbow_ramp(double X, double Lowval, double Upval);

/* glxwu_colour(): returns a const ptr to the RGB colour vector
 * indexed by Colidx. Colidx<0 returns black, Colidx>=GLXWU_COLNO
 * returns white.
 */
const GLfloat *glxwu_colour(int Colidx);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER glxwinutils.h ==== */
#endif	/* GLXWINUTILS_H */

#endif /* USE_OPENGL_GRAPHICS */
