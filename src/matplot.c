#ifdef USE_OPENGL_GRAPHICS

/* ==== PROJECT DRAGON: FUNCTIONS matplot.c ==== */

/* OpenGL program for visualising general rectangular MxN matrices. */

/* ANSI C+SGI OpenGL, 18-May-1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <GL/glu.h>

/* ---- MODULE HEADER ---- */

#include "matplot.h"

/* ==== FUNCTIONS ==== */

/* ---- Drawing ---- */

/* create_drawmat(): creates a Drawmatrix_ structure
 * which has a size RxC and the lower and upper limits
 * are set to Low and Up, respectively.
 * Returns a ptr to the newly allocated object.
 */
Drawmatrix_ *create_drawmat(unsigned int R, unsigned int C, 
	double Low, double Up)
{
    Drawmatrix_ *D=NULL;
    register unsigned int i;
    
    D=(Drawmatrix_*) malloc(sizeof(Drawmatrix_));
    D->Mat=(double **) calloc(R, sizeof(double *));  /* row ptrs */
    D->Mat[0]=(double *) calloc(R*C, sizeof(double));	/* elements */
    for (i=1; i<R; i++)
	D->Mat[i]=D->Mat[i-1]+C;    /* set other row ptrs */
    D->Row=R; D->Col=C;
    D->Lowval=Low; D->Upval=Up;
    D->Resizeneeded=GL_TRUE;	/* must calc scale */
    return(D);
}
/* END of create_drawmat() */

/* delete_drawmat(): frees up the storage associated with
 * *Drawmat, then removes the object itself.
 */
void delete_drawmat(Drawmatrix_* Drawmat)
{
    if (Drawmat!=NULL)
    {
	free(Drawmat->Mat[0]);
	free(Drawmat->Mat);
	free(Drawmat);
    }
}
/* END of delete_drawmat() */

/* plot_mat: creates a simple colour-coded dot representation of the
 * matrix stored in the Drawmat structure.
 * Plotting is performed in the window identified by *Windowinfo.
 */
void plot_mat(Windowinfo_ *Winfo, Drawmatrix_ *Drawmat)
{
    register unsigned int i,j;
    unsigned int Width=0, Height=0;
    XWindowAttributes Winattr;
    GLfloat Xscale=1.0, Yscale=1.0;

    /* connect the context to the window */
    glXMakeCurrent(Winfo->Dpy, Winfo->Win, Winfo->Ctx);
    
    /* Set the viewport to the full window
     * scale the matrix into it and then
     * project as 2D orthographic
     * if a resize is needed. Redrawing happens anyway
     */
    if (Drawmat->Resizeneeded)	/* first invocation */
    {
	XGetWindowAttributes(Winfo->Dpy, Winfo->Win, &Winattr);
	Width=Winattr.width;
	Height=Winattr.height;
    }
    else if (read_events(Winfo))    /* only if config changed */
    {
	if (Winfo->Event.type==ConfigureNotify)
	{
	    Width=Winfo->Event.xconfigure.width;
	    Height=Winfo->Event.xconfigure.height;
	    Drawmat->Resizeneeded=GL_TRUE;
	}
    }
    if (Drawmat->Resizeneeded)	/* modify projection and scaling */
    {
	Xscale=(GLfloat)Width/(GLfloat)(Drawmat->Col);
	Yscale=(GLfloat)Height/(GLfloat)(Drawmat->Row);
	
	glViewport(0, 0, Width, Height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, Width, 0.0, Height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(Xscale, Yscale, 1.0);
	Drawmat->Resizeneeded=GL_FALSE;
    }
    
    /* drawing into back buffer if double-buffer mode is on,
     * otherwise clear front buffer
     */
    if (Winfo->Dblbuffer)
	glDrawBuffer(GL_BACK);
    else    /* no double buffering :-( */
    {
	glDrawBuffer(GL_FRONT);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
    }
    
    for (i=0; i<Drawmat->Row; i++)           /* scans Mat */
	for (j=0; j<Drawmat->Col; j++)
	{
	    /* select code color for Mt values */
	    glColor3fv(rainbow_ramp(Drawmat->Mat[i][j], 
		Drawmat->Lowval, Drawmat->Upval));

	    /* plot the point */
	    glRecti(j, Drawmat->Row-i, j+1, Drawmat->Row-i-1);
 	}

    /* flush graphics */
    if (Winfo->Dblbuffer)
	glXSwapBuffers(Winfo->Dpy, Winfo->Win);
    else glFlush();
}
/* END of plot_mat */

/* ==== END OF FUNCTIONS matplot.c ==== */

#endif	/* USE_OPENGL_GRAPHICS */
