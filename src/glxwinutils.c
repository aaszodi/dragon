#ifdef USE_OPENGL_GRAPHICS

/* ==== FUNCTIONS glxwinutils.c ==== */

/* A few utility routines to open a GL window under X. */

/* ANSI C + X11 + OpenGL, 18-May-1998. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>

/* ---- MODULE HEADER ---- */

#include "glxwinutils.h"

/* ---- FILE_SCOPE VARS ---- */

/* Glxwucols: this array stores 3-component RGB colour
 * values which are indexed by Glxwucoltype_ enum constants.
 */
static const GLfloat Glxwucols[GLXWU_COLNO][3]={
		{0.0, 0.0, 0.0},    /* GLXWU_BLACK */
		{0.0, 0.0, 1.0},    /* GLXWU_BLUE */
		{0.0, 1.0, 1.0},    /* GLXWU_CYAN */
		{0.0, 1.0, 0.0},    /* GLXWU_GREEN */
		{1.0, 1.0, 0.0},    /* GLXWU_YELLOW */
		{1.0, 0.0, 0.0},    /* GLXWU_RED */
		{1.0, 1.0, 1.0}	    /* GLXWU_WHITE */
	    };

/*
 * Glxwudpy: we assume that all windows opened by the application
 * will be displayed on the same display. The pointer to
 * the corresponding X Display structure is stored here
 * after the first call to create_glxwindow().
 */
static Display *Glxwudpy=NULL;

/*
 * Openwinno: counts the number of open windows on *Glxwudpy.
 * After each call to create_glxwindow(), this is incremented by 1,
 * after destroy_glxwindow() is invoked, it is decremented by 1.
 * The display connection is destroyed if Openwinno==0.
 */
static int Openwinno=0;

/* EVENT_MASK: a constant for storing the flags of events
 * to be monitored.
 */
static const long EVENT_MASK = ExposureMask | StructureNotifyMask
	    	    	    | VisibilityChangeMask | KeyPressMask 
	    	    	    | ButtonPressMask | ButtonReleaseMask
	    	    	    | PointerMotionMask;

/* ---- PROTOTYPES ---- */

static Bool WaitForNotify(Display *d, XEvent *e, char *arg);

/* ==== FUNCTIONS ==== */

/* ---- Initialisation ---- */

/* create_glxwindow(): opens a GL/X window with its lower left
 * corner at Xorig:Yorig, size Xsize:Ysize, sets the title to Title and puts
 * all the gory details into *Windowinfo.
 * Aspect ratio is kept at Xsize/Ysize.
 * Returns GL_FALSE on error, GL_TRUE on success.
 */
GLenum create_glxwindow(Windowinfo_ *Winfo, int Xorig, int Yorig, 
	int Xsize, int Ysize, const char *Title)
{
    static int SnglAttrList[]={GLX_RGBA, 
	    GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None};
    static int DblAttrList[]={GLX_RGBA, GLX_DOUBLEBUFFER, 
	    GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None};
    XSizeHints *Shint=NULL;
    
    /* get local display connection: if already established, then just store */
    if (!Openwinno) 	/* no connection yet */
    {
    	Winfo->Dpy=Glxwudpy=XOpenDisplay(NULL);
	if (!Winfo->Dpy)
	{
	    fputs("\n! create_glxwindow(): Cannot open display\n", stderr);
	    return(GL_FALSE);
	}
	
	/* check if GLX is supported */
	if (!glXQueryExtension(Winfo->Dpy, NULL, NULL))
	{
	    fputs("\n! create_glxwindow(): GLX not supported\n", stderr);
	    return(GL_FALSE);
	}
    
	Openwinno=1;   /* OK, 1st window open */
    }
    else    /* connection already established */
    {
	Winfo->Dpy=Glxwudpy;
	Openwinno++;	    	/* count number of open windows */
    }
    
    /* get an appropriate visual: try double buffering first */
    Winfo->Dblbuffer=GL_TRUE;
    Winfo->Visinfo=glXChooseVisual(Winfo->Dpy, DefaultScreen(Winfo->Dpy), 
	    DblAttrList);
    if (Winfo->Visinfo==NULL)
    {
	fputs("\n? create_glxwindow(): Double buffering not supported\n", stderr);
	Winfo->Dblbuffer=GL_FALSE;
	
	/* trying single buffer */
	Winfo->Visinfo=glXChooseVisual(Winfo->Dpy, DefaultScreen(Winfo->Dpy), 
		SnglAttrList);
	if (Winfo->Visinfo==NULL)
	{
	    fputs("\n! create_glxwindow(): No suitable visual\n", stderr);
	    return(GL_FALSE);
	}
    }
    
    /* create GLX context */
    Winfo->Ctx=glXCreateContext(Winfo->Dpy, Winfo->Visinfo, 0, GL_TRUE);
    if (Winfo->Ctx==NULL)
    {
	fputs("\n! create_glxwindow(): Cannot create GL context\n", stderr);
	return(GL_FALSE);
    }
    
    /* create colormap */
    Winfo->Winattr.colormap=XCreateColormap(Winfo->Dpy,
	    RootWindow(Winfo->Dpy, Winfo->Visinfo->screen), 
	    Winfo->Visinfo->visual, AllocNone);
    
    /* create the window */
    Winfo->Winattr.border_pixel=0;
    Winfo->Winattr.event_mask=EVENT_MASK;
    Winfo->Win=XCreateWindow(Winfo->Dpy, RootWindow(Winfo->Dpy, Winfo->Visinfo->screen), 
	    Xorig, Yorig, Xsize, Ysize, 1, 
	    Winfo->Visinfo->depth, InputOutput, Winfo->Visinfo->visual, 
	    CWBorderPixel|CWColormap|CWEventMask, &(Winfo->Winattr));
    
    /* add the title */
    XStoreName(Winfo->Dpy, Winfo->Win, Title);
    
    /* do the size hints */
    Shint=XAllocSizeHints();
    Shint->flags=USPosition|USSize|PMinSize|PAspect;
    Shint->min_width=100;   /* minimal sizes */
    Shint->min_height=100;
    Shint->min_aspect.x=Shint->max_aspect.x=Xsize;
    Shint->min_aspect.y=Shint->max_aspect.y=Ysize;
    XSetWMProperties(Winfo->Dpy, Winfo->Win, NULL, NULL, 
	    NULL, 0, Shint, NULL, NULL);
    XFree(Shint);
    
    XMapWindow(Winfo->Dpy, Winfo->Win);
    XIfEvent(Winfo->Dpy, &(Winfo->Event), WaitForNotify, (char *)(Winfo->Win));
    
    /* connect the context to the window */
    glXMakeCurrent(Winfo->Dpy, Winfo->Win, Winfo->Ctx);
    return(GL_TRUE);    
}
/* END of create_glxwindow() */

/* WaitForNotify: needed by XIfEvent. */
static Bool WaitForNotify(Display *d, XEvent *e, char *arg)
{
     return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}

/* destroy_glxwindow(): removes the window *Winfo from the screen. */
void destroy_glxwindow(Windowinfo_ *Winfo)
{
    if (Winfo==NULL || Winfo->Dpy==NULL) return;
    
    glFlush(); glFinish();
    glXMakeCurrent(Winfo->Dpy,None,NULL);   /* no current context */
    glXDestroyContext(Winfo->Dpy, Winfo->Ctx);
    XFreeColormap(Winfo->Dpy, Winfo->Winattr.colormap);
    XFree((char *)(Winfo->Visinfo));
    XDestroyWindow(Winfo->Dpy, Winfo->Win);
    Winfo->Dpy = NULL;
    Openwinno--;    	/* one window less */
    if (!Openwinno) 	/* break X server connection as well */
    {
	XCloseDisplay(Glxwudpy);
	Glxwudpy=NULL;
    }
}
/* END of destroy_glxwindow() */

/* ---- Event handling ---- */

/* read_events(): reads the event queue associated with Winfo->Dpy "dry"
 * and returns GL_TRUE if Expose, VisibilityNotify or ConfigureNotify
 * events were detected. If only Expose was seen then Winfo->Event.type will
 * contain Expose, if at least one ConfigureNotify was seen then
 * Winfo->Event will contain the full description of it (to help
 * with resizing).
 */
GLenum read_events(Windowinfo_ *Winfo)
{
    XEvent Event;
    GLenum Redraw=GL_FALSE;
    
    Event.type=Winfo->Event.type=0;
    /* collect X events and read the queue dry
     * while condition was XPending(Winfo->Dpy)
     */
    while(XCheckWindowEvent(Winfo->Dpy, Winfo->Win, EVENT_MASK, &Event))
    {
	switch (Event.type)
	{
	    case Expose:
	    case VisibilityNotify:
		if (Winfo->Event.type!=ConfigureNotify)
		{
		    Winfo->Event=Event;	/* save only if no ConfigureNotify was seen */
		    Redraw=GL_TRUE;
		}
	    break;
	    case ConfigureNotify:
		Winfo->Event=Event; /* overrides whatever was there before */
		Redraw=GL_TRUE;
	    break;
	    default: break;
	}	/* switch */
    }
    return(Redraw);
}
/* END of read_events() */

/* ---- Colour ---- */

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
const GLfloat *rainbow_ramp(double X, double Lowval, double Upval)
{
    static GLfloat Color[3];

    if (Upval==Lowval)
    {
	Color[0]=Color[1]=Color[2]=0.0;	    /* pitch black */
	return((const GLfloat*)Color);
    }
    
    /* rescale: Lowval=0, Upval=1 */
    X=(X-Lowval)/(Upval-Lowval);
    
    /* find ranges, set colours */
    if (X<0.0)	/* out of range */
    {
	Color[0]=Color[1]=Color[2]=0.0;	    /* pitch black */
    }    
    else if (X<=0.25)    /* blue...cyan: full blue+green grows */
    {
	Color[0]=0.0;	/* no red */
	Color[1]=4.0*X;	/* green grows */
	Color[2]=1.0;	/* full blue */
    }
    else if (X<=0.5)	/* cyan..green: full green+blue shrinks */
    {
	Color[0]=0.0;	/* no red */
	Color[1]=1.0;	/* full green */
	Color[2]=4.0*(0.5-X);	/* blue shrinks */
    }
    else if (X<=0.75)
    {
	Color[0]=4.0*(X-0.5);	/* red grows */
	Color[1]=1.0;	/* full green */
	Color[2]=0.0;	/* no blue */
    }
    else if (X<=1.0)
    {
	Color[0]=1.0;	/* full red */
	Color[1]=4.0*(1.0-X);	/* green shrinks */
	Color[2]=0.0;	/* no blue */
    }
    else    /* out of range */
    {
	Color[0]=Color[1]=Color[2]=1.0;	    /* snow white */
    }
    return((const GLfloat*)Color);
}
/* END of rainbow_ramp */

/* glxwu_colour(): returns a const ptr to the RGB colour vector
 * indexed by Colidx. Colidx<0 returns black, Colidx>=GLXWU_COLNO
 * returns white.
 */
const GLfloat *glxwu_colour(int Colidx)
{
    if (Colidx<0) Colidx=0;
    if (Colidx>=GLXWU_COLNO) Colidx=GLXWU_COLNO-1;
    return((const GLfloat*)Glxwucols[Colidx]);
}
/* END of glxwu_colour() */

/* ==== END OF FUNCTIONS glxwinutils.c ==== */

#endif	/* USE_OPENGL_GRAPHICS */
