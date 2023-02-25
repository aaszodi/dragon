# Makefile for DRAGON-IV
# 27-Aug-1998. (c) Andras Aszodi

# The directory layout is as follows:-
#
# dragon4: root to the whole project
# |
# +-- src++: *.c++, *.h files for the C++ modules
# |
# +-- src: *.[ch] files for the C modules
# |
# +-- test: c++, *.[ch] files for the module test programs
# |
# +-- bin: objects,libraries,executables <--YOU ARE HERE
# |    |
# |    +-- irix-o32: old 32-bit SGI ABI     \
# |    |                                    |
# |    +-- irix-n32: new 32/64-bit SGI ABI  |
# |    |                                    |
# |    +-- irix-n64: new 64-bit SGI ABI     |
# |    |                                    |-- $(ABI) directories
# |    +-- solaris-sparc: GCC/Solaris SPARC | 
# |    |                                    |
# |    +-- solaris-i86 : GCC/Solaris x86    | 
# |    |                                    |
# |    +-- linux: GCC/Linux ABI             /
# |
# +-> lib: root of the utility hierarchy (symlink to ~/lib)
#      |
#      +-- irix-o32: old 32-bit SGI ABI     \
#      |                                    |
#      +-- irix-n32: new 32/64-bit SGI ABI  |
#      |                                    |
#      +-- irix-n64: new 64-bit SGI ABI     |
#      |                                    |-- $(UTILS) directories
#      +-- solaris-sparc: GCC/solaris SPARC |
#      |                                    |
#      +-- solaris-i86 : GCC/Solaris x86    | 
#      |                                    |
#      +-- linux: GCC/Linux ABI             /
#      |
#      +-- c: C source tree
#      |   |
#      |   +-- incl: header files (*.h)
#      |   |
#      |   +-- src: function definition files (*.c)
#      |
#      +-- cc: C++ source tree
#          |
#          +-- incl: non-template class header files (*.h)
#          |
#          +-- src: non-template function definition files (*.c++)
#          |
#          +-- tmpl: template declarations (*.h) and definitions (*.c++)

# ---- ABI-DEPENDENT DEFINITIONS ----

# NOTE: this Makefile supports the maintenance of incompatible
# ABIs. Currently SGI o32,n32,n64,GCC/Linux,GCC/Solaris are supported.
# The Top Makefile should cd to bin/$(ABI) and invoke make from
# there as make -f ../Makefile (which is this file).
# ABI: currently (irix-(o32|n32|n64)|solaris-(sparc|i86)|linux),
# selects the application binary interface directory

# ---- Local settings ----

# Note: macros prefixed with "CC" are for the C++ compiler,
# the macros beginning with "C" are for the C compiler.

# Source directories
CSRC = ../../src
CCSRC = ../../src++

# Library
LIB = ../../lib
UTILS = $(LIB)/$(ABI)

# Utility header locations
CHDR = $(LIB)/c/incl
CCHDR = $(LIB)/cc/incl
TMPLHDR = $(LIB)/cc/tmpl

# ---- ABI-specific macros ----

# The settings for a given ABI or architecture come from two
# include files which are called "Makefile.<ABI>" in the "lib" 
# and "bin" directories. The first contains general settings
# for the architecture and is supposed to be used by all programs
# using the utilities in "lib". The second contains settings
# specific to the DRAGON package (linking etc.).
# When porting to a new ABI, new include files must be prepared.

include $(LIB)/Makefile.$(ABI)
include ../Makefile.$(ABI)

PVMHDR = $(PVMROOT)/include

# ---- All compiler flags together ----

CFLAGS = $(LIBABIFG) $(BINABIFG) $(COPT) -I$(CHDR)
CCFLAGS = $(LIBABIFG) $(BINABIFG) $(CCOPT) -I$(CCHDR) -I$(TMPLHDR) -I$(PVMHDR)

# ==== C++ MODULES ====

# ---- LIBRARIES ----

# ---- libpoly.a ----

# the model chain class (alignment, geometry etc.) collection
libpoly.a: libpoly.a(Acdist.o) libpoly.a(Align.o) libpoly.a(Distpred.o) \
		libpoly.a(Polymer.o) libpoly.a(Property.o) libpoly.a(Simil.o)
	rm $?

# Side-chain atom distances from centroid and C-alpha
libpoly.a(Acdist.o): $(CCSRC)/Acdist.c++ $(CCSRC)/Acdist.h $(TMPLHDR)/Array.h $(CCHDR)/String.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Acdist.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Multiple alignment storage
libpoly.a(Align.o): $(CCSRC)/Align.c++ $(CCSRC)/Align.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Align.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Distance prediction
libpoly.a(Distpred.o): $(CCSRC)/Distpred.c++ $(CCSRC)/Distpred.h \
		$(TMPLHDR)/Array.h $(CCHDR)/Vector.h $(CCHDR)/Spl.h $(CCHDR)/Cdf.h $(CCHDR)/Pmest.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Distpred.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Polymer data
libpoly.a(Polymer.o): $(CCSRC)/Polymer.c++ $(CCSRC)/Polymer.h \
		$(CCSRC)/Acdist.h $(CCSRC)/Align.h $(CCSRC)/Distpred.h \
		$(CCSRC)/Property.h $(CCSRC)/Simil.h \
		$(CCHDR)/Stat2.h $(TMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Polymer.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Monomer properties
libpoly.a(Property.o): $(CCSRC)/Property.c++ $(CCSRC)/Property.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Property.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Similarity matrix
libpoly.a(Simil.o): $(CCSRC)/Simil.c++ $(CCSRC)/Simil.h $(CCHDR)/Sqmat.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Simil.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# ---- libpieces.a ----

# the secondary-structure segment handling class collection

libpieces.a: libpieces.a(Beta.o) libpieces.a(Helix.o) libpieces.a(Pieces.o) \
		libpieces.a(Segment.o) libpieces.a(Sstrbase.o) 
	rm $?

# beta-sheets
libpieces.a(Beta.o): $(CCSRC)/Beta.c++ $(CCSRC)/Beta.h \
		$(CCSRC)/Sstrbase.h $(CCSRC)/Segment.h \
		$(CCHDR)/Bits.h $(CCHDR)/Points.h $(TMPLHDR)/Array.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Beta.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# alpha-helices
libpieces.a(Helix.o): $(CCSRC)/Helix.c++ $(CCSRC)/Helix.h \
		$(CCSRC)/Sstrbase.h $(CCSRC)/Segment.h \
		$(CCHDR)/Bits.h $(CCHDR)/Points.h $(TMPLHDR)/Array.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Helix.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Pieces of secstr together
libpieces.a(Pieces.o): $(CCSRC)/Pieces.c++ $(CCSRC)/Pieces.h \
		$(CCSRC)/Segment.h $(CCSRC)/Sstrbase.h \
		$(CCSRC)/Helix.h $(CCSRC)/Beta.h $(CCHDR)/Bits.h \
		$(TMPLHDR)/Array.h $(TMPLHDR)/List1.h $(TMPLHDR)/Smartptr.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Pieces.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Chain segments
libpieces.a(Segment.o): $(CCSRC)/Segment.c++ $(CCSRC)/Segment.h \
		$(CCHDR)/Bits.h $(TMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Segment.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Secondary structure base class
libpieces.a(Sstrbase.o): $(CCSRC)/Sstrbase.c++ $(CCSRC)/Sstrbase.h $(CCSRC)/Segment.h \
		$(CCHDR)/Bits.h $(CCHDR)/Points.h $(TMPLHDR)/Array.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Sstrbase.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# ---- libgraph.a ----

# OpenGL graphics
libgraph.a: libgraph.a(Graphics.o) libgraph.a(cadraw.o) \
    	libgraph.a(matplot.o) libgraph.a(glxwinutils.o)
	rm $?

# Graphics class
libgraph.a(Graphics.o): $(CCSRC)/Graphics.c++ $(CCSRC)/Graphics.h $(CCSRC)/Polymer.h \
		$(CCHDR)/Points.h $(CSRC)/cadraw.h $(CSRC)/matplot.h $(CSRC)/glxwinutils.h
	$(CXX) $(CCFLAGS) -I$(CSRC) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/Graphics.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Molecule draw
libgraph.a(cadraw.o): $(CSRC)/cadraw.c $(CSRC)/cadraw.h $(CSRC)/glxwinutils.h
	$(CC) $(CFLAGS) -c $(CSRC)/cadraw.c -o $%
	$(AR) $(ARFLAGS) $@ $%

# Distance matrix plot
libgraph.a(matplot.o): $(CSRC)/matplot.c $(CSRC)/matplot.h $(CSRC)/glxwinutils.h
	$(CC) $(CFLAGS) -c $(CSRC)/matplot.c -o $%
	$(AR) $(ARFLAGS) $@ $%

# OpenGL/X window utilities
libgraph.a(glxwinutils.o): $(CSRC)/glxwinutils.c $(CSRC)/glxwinutils.h
	$(CC) $(CFLAGS) -c $(CSRC)/glxwinutils.c -o $%
	$(AR) $(ARFLAGS) $@ $%

# ---- MAIN MODULES ----

# DRAGON
Dragon.o: $(CCSRC)/Dragon.c++ $(CCSRC)/Access.h $(CCSRC)/Clip.h $(CCSRC)/Density.h \
		$(CCSRC)/Graphics.h $(CCSRC)/Hmom.h $(CCSRC)/Homodel.h $(CCSRC)/Iproj.h \
		$(CCSRC)/Params.h $(CCSRC)/Pieces.h $(CCSRC)/Polymer.h $(CCSRC)/Pvmtask.h $(CCSRC)/Output.h \
		$(CCSRC)/Restr.h $(CCSRC)/Score.h $(CCSRC)/Sigproc.h $(CCSRC)/Steric.h \
		$(CCSRC)/Sterchem.h $(CCSRC)/Tangles.h $(CCSRC)/Viol.h \
		$(CSRC)/version.h $(CHDR)/cmdopt.h $(CHDR)/tstamp.h $(CCHDR)/String.h
	$(CXX) $(CCFLAGS) -I$(CSRC) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/Dragon.c++ -o $@ 

# Quanta H-bond assignment
hbassign.o: $(CCSRC)/hbassign.c++ $(CCSRC)/Align.h $(CHDR)/dsspread.h $(TMPLHDR)/List1.h
	$(CXX) $(CCFLAGS) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/hbassign.c++ -o $@

# Secondary structure mapping
Secmap.o: $(CCSRC)/Secmap.h $(CCSRC)/Secmap.c++
	$(CXX) $(CCFLAGS) -I$(CHDR) -c $(CCSRC)/Secmap.c++ -o $@
secmap.o: $(CCSRC)/secmap.c++ $(CCSRC)/Align.h $(CCSRC)/Secmap.h $(CHDR)/dsspread.h
	$(CXX) $(CCFLAGS) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/secmap.c++ -o $@

# side chain homology
Aacid.o: $(CCSRC)/Aacid.c++ $(CCSRC)/Aacid.h $(CCHDR)/Points.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Aacid.c++ -o $@
sidech.o: $(CCSRC)/sidech.c++ $(CCSRC)/Align.h $(CCSRC)/Aacid.h \
		$(CCHDR)/Hirot.h $(TMPLHDR)/Array.h $(CHDR)/pdbprot.h
	$(CXX) $(CCFLAGS) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/sidech.c++ -o $@

# ---- MODULES ----

# Accessibility
Access.o: $(CCSRC)/Access.c++ $(CCSRC)/Access.h $(CCSRC)/Fakebeta.h $(CCSRC)/Polymer.h \
		$(CCHDR)/Bits.h $(CCHDR)/Points.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Access.c++ -o $@

# Command line interpreter
Clip.o: $(CCSRC)/Clip.c++ $(CCSRC)/Clip.h $(CCSRC)/Params.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Clip.c++ -o $@

# Density
Density.o: $(CCSRC)/Density.c++ $(CCSRC)/Density.h $(CCSRC)/Pieces.h \
		$(CCHDR)/Points.h $(CCHDR)/Trimat.h \
		$(CCHDR)/Sqmat.h $(CCHDR)/Ql.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Density.c++ -o $@

# Fake CB distances
Fakebeta.o: $(CCSRC)/Fakebeta.c++ $(CCSRC)/Fakebeta.h $(CCSRC)/Polymer.h $(CCHDR)/Points.h \
		$(TMPLHDR)/Array.h $(CCHDR)/Sqmat.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Fakebeta.c++ -o $@

# Hydrophobic moments
Hmom.o: $(CCSRC)/Hmom.c++ $(CCSRC)/Hmom.h $(CCSRC)/Fakebeta.h $(CCSRC)/Pieces.h $(CCSRC)/Polymer.h \
		$(CCHDR)/Points.h $(CCHDR)/Svd.h \
		 $(CCHDR)/Vector.h $(CCHDR)/Sqmat.h 
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Hmom.c++ -o $@

# Homology modelling
Homodel.o: $(CCSRC)/Homodel.c++ $(CCSRC)/Homodel.h $(CCSRC)/Polymer.h \
		$(CCSRC)/Pvmtask.h $(CCSRC)/Restr.h \
		$(CCHDR)/Points.h $(CCHDR)/Hirot.h $(CCHDR)/Vector.h $(CHDR)/pdbprot.h
	$(CXX) $(CCFLAGS) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/Homodel.c++ -o $@

# Inertial projection
Iproj.o: $(CCSRC)/Iproj.c++ $(CCSRC)/Iproj.h $(CCSRC)/Pieces.h $(CCHDR)/Rsmdiag.h \
		$(CCHDR)/Vector.h $(CCHDR)/Trimat.h $(CCHDR)/Sqmat.h \
		$(CCHDR)/Hirot.h $(CCHDR)/Ql.h $(CCHDR)/Points.h $(TMPLHDR)/Array.h $(TMPLHDR)/Maskarr.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Iproj.c++ -o $@

# Output to PDB 
Output.o: $(CCSRC)/Output.c++ $(CCSRC)/Output.h $(CCSRC)/Fakebeta.h \
		$(CCSRC)/Polymer.h $(CCSRC)/Pieces.h $(CCSRC)/Score.h $(CSRC)/version.h \
		$(CCHDR)/String.h $(TMPLHDR)/List1.h $(CHDR)/pdbprot.h
	$(CXX) $(CCFLAGS) -I$(CSRC) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/Output.c++ -o $@

# Parameters
Paramstr.o: $(CCSRC)/Paramstr.c++ $(CCSRC)/Paramstr.h $(CCSRC)/Parambase.h \
		$(CCHDR)/String.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Paramstr.c++ -o $@

Params.o: $(CCSRC)/Params.c++ $(CCSRC)/Params.h \
		$(CCSRC)/Parambase.h $(CCSRC)/Paramstr.h $(CCSRC)/Paramlim.h \
		$(CCHDR)/String.h $(TMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Params.c++ -o $@

# Parallel Virtual Machine support
Pvmtask.o: $(CCSRC)/Pvmtask.c++ $(CCSRC)/Pvmtask.h $(CCSRC)/Params.h \
		$(CCSRC)/Sigproc.h $(CCHDR)/String.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Pvmtask.c++ -o $@

# Distance restraints
Restr.o: $(CCSRC)/Restr.c++ $(CCSRC)/Restr.h \
		$(CCSRC)/Polymer.h $(CCSRC)/Pieces.h $(CHDR)/portrandom.h $(CCHDR)/Bits.h $(TMPLHDR)/Array.h $(TMPLHDR)/List1.h
	$(CXX) $(CCFLAGS) -I$(CHDR) $(TMPLOPTS) -c $(CCSRC)/Restr.c++ -o $@

# Scores
Score.o: $(CCSRC)/Score.c++ $(CCSRC)/Score.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Score.c++ -o $@

# Signal traps and multiple process management
Sigproc.o: $(CCSRC)/Sigproc.c++ $(CCSRC)/Sigproc.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Sigproc.c++ -o $@

# Spectral Gradient
Specgrad.o: $(CCSRC)/Specgrad.c++ $(CCSRC)/Specgrad.h $(CCHDR)/Points.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Specgrad.c++ -o $@

# Steric adjustments
Steric.o: $(CCSRC)/Steric.c++ $(CCSRC)/Steric.h \
		$(CCSRC)/Pieces.h $(CCSRC)/Restr.h $(CCSRC)/Polymer.h \
		$(CCSRC)/Fakebeta.h $(CCSRC)/Specgrad.h \
		$(CCSRC)/Score.h $(CCSRC)/Viol.h \
	$(CCHDR)/Trimat.h $(CCHDR)/Vector.h $(CCHDR)/Points.h $(CCHDR)/Hirot.h $(TMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Steric.c++ -o $@

# General stereochemical adjustments
Sterchem.o: $(CCSRC)/Sterchem.c++ $(CCSRC)/Sterchem.h $(CCSRC)/Pieces.h \
		$(CCHDR)/Points.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Sterchem.c++ -o $@

# Tangle detection and elimination
Tangles.o: $(CCSRC)/Tangles.c++ $(CCSRC)/Tangles.h $(CCSRC)/Pieces.h \
		$(CCHDR)/Bits.h $(CCHDR)/Points.h $(CCHDR)/Svd.h \
		$(TMPLHDR)/Array.h $(TMPLHDR)/List1.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Tangles.c++ -o $@

# Violation list
Viol.o: $(CCSRC)/Viol.c++ $(CCSRC)/Viol.h $(TMPLHDR)/List1.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Viol.c++ -o $@

# ==== C MODULES ====

# ---- rank ----

# ranks the output files according to score
rank.o: $(CSRC)/rank.c $(CHDR)/cmdopt.h
	$(CC) $(CFLAGS) -c $(CSRC)/rank.c -o $@

# ---- clumsy ----

# Similarity with hierarchic metric single-linkage clustering
clumsy.o: $(CSRC)/clumsy.c $(CSRC)/dslclu.h $(CSRC)/rotpdb.h \
		 $(CHDR)/cmdopt.h $(CHDR)/siva.h
	$(CC) $(CFLAGS) -c $(CSRC)/clumsy.c -o $@

# PDB I/O and optimal rotation
rotpdb.o: $(CSRC)/rotpdb.c $(CSRC)/rotpdb.h \
		$(CHDR)/matrix.h $(CHDR)/bestrot.h $(CHDR)/pdbprot.h 
	$(CC) $(CFLAGS) -c $(CSRC)/rotpdb.c -o $@

# single-linked metric distance clustering
dslclu.o: $(CSRC)/dslclu.c $(CSRC)/dslclu.h $(CHDR)/matrix.h
	$(CC) $(CFLAGS) -c $(CSRC)/dslclu.c -o $@

# ==== PROGRAMS ====

# C++ utility objects
C++UTILOBJS =  $(UTILS)/Hirot.o $(UTILS)/Points.o

# C utility objects
CUTILOBJS = $(UTILS)/cmdopt.o $(UTILS)/pdbprot.o \
	$(UTILS)/portrandom.o $(UTILS)/tstamp.o

# C++ modules
C++MODOBJS = Access.o Clip.o Density.o \
		Fakebeta.o Hmom.o Homodel.o \
		Iproj.o Output.o Paramstr.o Params.o Pvmtask.o \
		Restr.o Score.o Sigproc.o \
		Specgrad.o Steric.o \
		Sterchem.o Tangles.o Viol.o

# Library paths
LIBPATH = -L$(UTILS) -L../$(ABI) $(PVMPATH)

# Libraries
LIBS = -lpoly -lpieces -lccstat -lccutils -linalg  $(PVMLIB)  $(GLIB) $(MATHLIB)

# DRAGON
dragon: $(UTILS)/libinalg.a $(UTILS)/libccstat.a $(UTILS)/libccutils.a $(C++UTILOBJS) \
		libpoly.a libpieces.a libgraph.a \
		Dragon.o $(C++MODOBJS) $(C++UTILOBJS) $(CUTILOBJS)
	$(CC) $(CFLAGS) -DDRAGON_VERSION='"$(VERSION)-$(ABI)"' -c $(CSRC)/version.c -o version.o;
	$(CXX) $(CCFLAGS) $(TMPLOPTS) $(C++UTILOBJS) Dragon.o $(C++MODOBJS) \
		$(CUTILOBJS) version.o \
		$(LIBPATH) $(LIBS) -o $@
# Use this for Linux memory leak detection after $(LIBS)
# /usr/local/lib/ccmalloc.o -ldl
#
# Use this for SGI memory leak detection (if you have CaseVision/Workshop)
# -lmalloc_cv
#

# QUANTA H-bond restraint generation
hbassign: hbassign.o libpoly.a $(UTILS)/dsspread.o
	$(CXX) $(CCFLAGS) $(TMPLOPTS) hbassign.o $(UTILS)/dsspread.o \
		-L../$(ABI) -lpoly -lm -o $@

# Secondary structure mapping
secmap: secmap.o Secmap.o libpoly.a $(UTILS)/dsspread.o
	$(CXX) $(CCFLAGS) $(TMPLOPTS) secmap.o Secmap.o $(UTILS)/dsspread.o \
		-L../$(ABI) -lpoly -lm -o $@

# Homology-based side chains
sidech: sidech.o Aacid.o $(UTILS)/Hirot.o $(UTILS)/pdbprot.o \
		$(UTILS)/libccutils.a $(UTILS)/libinalg.a
	$(CXX) $(CCFLAGS) $(TMPLOPTS) sidech.o \
		Aacid.o $(UTILS)/Hirot.o $(UTILS)/Points.o $(UTILS)/pdbprot.o \
		-L$(UTILS) -L../$(ABI) -lpoly -lccutils -linalg -lm -o $@

# Output file ranking (C only)
rank: $(UTILS)/cmdopt.o rank.o
	$(CC) $(CFLAGS) $(UTILS)/cmdopt.o rank.o -o $@

# Structure clustering (C only)
clumsy: clumsy.o dslclu.o rotpdb.o \
		$(UTILS)/cmdopt.o $(UTILS)/bestrot.o $(UTILS)/pdbprot.o \
		$(UTILS)/siva.o $(UTILS)/matrix.o $(UTILS)/ql.o
	$(CC) $(CFLAGS) clumsy.o dslclu.o rotpdb.o $(UTILS)/cmdopt.o $(UTILS)/bestrot.o \
		$(UTILS)/pdbprot.o $(UTILS)/siva.o $(UTILS)/matrix.o $(UTILS)/ql.o \
		-lm -o $@
