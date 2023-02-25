# Description of the "lib" directory
# ==================================
# 
# This directory contains all my utilities, both C and C++.
# Since the introduction of the IRIX 6.2/MIPSPro 7.0 compilers,
# three application binary interfaces (ABIs) are supported on
# SGI machines. O32 is the old 32-bit ABI (MIPS I,II instruction set),
# N32 is the new 32-bit ABI with more floating-point registers,
# MIPS III/IV instruction set and 64-bit calling conventions;
# N64 is the fully 64-bit ABI (MIPS III/IV on R8000 and R10000 only).
# Object files and libraries for the different ABIs are stored
# separately as they cannot be linked together.
# In addition, the GCC/Linux, GCC/Solaris ABI is supported now.
# 
# Directory layout
# ================
#
# lib: root of the utility hierarchy  <--YOU ARE HERE
#  |
#  +-- irix-o32: old 32-bit SGI ABI     \
#  |                                    |
#  +-- irix-n32: new 32/64-bit SGI ABI  |
#  |                                    |
#  +-- irix-n64: new 64-bit SGI ABI     |
#  |                                    |-- $(UTILS) directories
#  +-- solaris-sparc: GCC/solaris SPARC |   Compilation performed here
#  |                                    |
#  +-- solaris-i86 : GCC/Solaris x86    | 
#  |                                    |
#  +-- linux: GCC/Linux ABI             /
#  |
#  +-- c: C source tree
#  |   |
#  |   +-- incl: header files (*.h)
#  |   |
#  |   +-- src: function definition files (*.c)
#  |
#  +-- cc: C++ source tree
#      |
#      +-- incl: non-template class header files (*.h)
#      |
#      +-- src: non-template function definition files (*.c++)
#      |
#      +-- tmpl: template declarations (*.h) and definitions (*.c++)

# Usage:
# The library modules are usually made by Top Makefiles in the
# project directories. Upon invocation, these should specify
# the target list, cd to the corresponding lib/$(ABI) subdirectory
# and invoke this Makefile from there as "make -f ../Makefile ABI='$(ABI)' ..."
# The following macro is expected from the Top Makefile:
# ABI: currently (irix-(o32|n32|n64)|solaris-(sparc|i86)|linux),
# selects the application binary interface directory
# As a convention, C++ macro names begin with "CC", C macros with "C".

# The shell
SHELL = /bin/sh

# ---- Subdirectories ----

# Source trees
CTREE = ../c
CCTREE = ../cc

# Header directories
CHDR = $(CTREE)/incl
CCHDR = $(CCTREE)/incl
CCTMPLHDR = $(CCTREE)/tmpl

# Source directories
CSRC = $(CTREE)/src
CCSRC = $(CCTREE)/src

# ---- ABI-specific macros ----

# The settings for a given ABI or architecture come from
# an include file which is called "Makefile.<ABI>" in the "lib" directory.
# When porting to a new ABI, a new include file must be prepared.

include ../Makefile.$(ABI)

# ---- All compiler flags together ----

CFLAGS = $(COPT) $(LIBABIFG) -I$(CHDR)
CCFLAGS = $(CCOPT) $(LIBABIFG) -I$(CCHDR) -I$(CCTMPLHDR) -I$(CHDR)

# ---- TARGET LIST ---- 

# When adding a new object or library to this Makefile,
# don't forget to update these macros below!

COBJECTS = bestrot.o \
	cmdopt.o \
	dsspread.o \
	matrix.o \
	portrandom.o \
	pdbprot.o \
	ql.o \
	siva.o \
	tstamp.o

CCLIBRARIES = libinalg.a libccstat.a libccutils.a
CCDSOS = $(LIBRARIES:.a=.so)
CCOBJECTS = Hirot.o Points.o

# ---- MAIN RULES ----

# generate all objects and libraries (default)
# add 'dsos' if DSOs are needed
all-c: objs-c
all-c-nosgi: objs-c
all-cc: libs-cc objs-cc

# make libraries only
libs-cc: $(CCLIBRARIES)

# make DSOs
# dsos-cc: $(CCDSOS)

# make object modules only
objs-c: $(COBJECTS)
objs-cc: $(CCOBJECTS)

# cleanup
clean-c:
	rm -f $(COBJECTS)
clean-cc:
	rm -f $(CCLIBRARIES) $(CCDSOS) $(CCOBJECTS)

# ---- NON-LIBRARY C MODULES ----

# McLachlan's best rotation
bestrot.o: $(CSRC)/bestrot.c $(CHDR)/bestrot.h $(CHDR)/matrix.h $(CHDR)/siva.h
	$(CC) $(CFLAGS) -c $(CSRC)/bestrot.c -o $@

# command-line options
cmdopt.o: $(CSRC)/cmdopt.c $(CHDR)/cmdopt.h
	$(CC) $(CFLAGS) -c $(CSRC)/cmdopt.c -o $@

# DSSP file processing
dsspread.o: $(CSRC)/dsspread.c $(CHDR)/dsspread.h
	$(CC) $(CFLAGS) -c $(CSRC)/dsspread.c -o $@

# Simple matrices
matrix.o: $(CSRC)/matrix.c $(CHDR)/matrix.h
	$(CC) $(CFLAGS) -c $(CSRC)/matrix.c -o $@

# Portable random number generators
portrandom.o: $(CSRC)/portrandom.c $(CHDR)/portrandom.h
	$(CC) $(CFLAGS) -c $(CSRC)/portrandom.c -o $@

# Protein-oriented PDB I/O
pdbprot.o: $(CSRC)/pdbprot.c $(CHDR)/pdbprot.h
	$(CC) $(CFLAGS) -c $(CSRC)/pdbprot.c -o $@

# QL decomposition
ql.o: $(CSRC)/ql.c $(CHDR)/ql.h $(CHDR)/matrix.h
	$(CC) $(CFLAGS) -c $(CSRC)/ql.c -o $@

# Singular Value Decomposition
siva.o: $(CSRC)/siva.c $(CHDR)/siva.h $(CHDR)/matrix.h $(CHDR)/ql.h
	$(CC) $(CFLAGS) -c $(CSRC)/siva.c -o $@

# Time stamps
tstamp.o: $(CSRC)/tstamp.c $(CHDR)/tstamp.h
	$(CC) $(CFLAGS) -c $(CSRC)/tstamp.c -o $@

# ---- C++ LIBRARIES ----

# Both static libraries (*.a) and DSOs (*.so) are provided.
# DSO building can be done via the "CC" driver (as opposed to C-language
# DSOs which are built using "ld").

# ---- libinalg.a ----

# The Linear Algebra library
libinalg.a: libinalg.a(Matbase.o) libinalg.a(Rectbase.o) \
		libinalg.a(Sqbase.o) libinalg.a(Matrix.o) \
		libinalg.a(Sqmat.o) libinalg.a(Trimat.o) \
		libinalg.a(Lu.o) libinalg.a(Ql.o) libinalg.a(Rsmdiag.o) \
		libinalg.a(Safety.o) libinalg.a(Svd.o) \
		libinalg.a(Vector.o) libinalg.a(Vmutils.o)
	rm $?

# DSO version
# libinalg.so: libinalg.a
#	cd $(ABI); $(CXX) $(SOFLAGS) -shared -rdata_shared -all $(@F:.so=.a) -o $(@F)

# ABC of all matrices
libinalg.a(Matbase.o): $(CCSRC)/Matbase.c++ $(CCHDR)/Matbase.h $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Matbase.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# ABC of rectangular matrices (Matrix,Sqmat)
libinalg.a(Rectbase.o): $(CCSRC)/Rectbase.c++ $(CCHDR)/Rectbase.h $(CCHDR)/Matbase.h $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Rectbase.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# ABC of square matrices (Sqmat,Trimat)
libinalg.a(Sqbase.o): $(CCSRC)/Sqbase.c++ $(CCHDR)/Sqbase.h $(CCHDR)/Matbase.h $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Sqbase.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# General RxC matrices
libinalg.a(Matrix.o): $(CCSRC)/Matrix.c++ $(CCHDR)/Matrix.h $(CCHDR)/Matbase.h $(CCHDR)/Rectbase.h $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Matrix.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# RxR square matrices
libinalg.a(Sqmat.o): $(CCSRC)/Sqmat.c++ $(CCHDR)/Sqmat.h \
		$(CCHDR)/Matbase.h $(CCHDR)/Rectbase.h $(CCHDR)/Sqbase.h $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Sqmat.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# symmetric RxR matrices stored in triangular form
libinalg.a(Trimat.o): $(CCSRC)/Trimat.c++ $(CCHDR)/Trimat.h \
		$(CCHDR)/Matbase.h $(CCHDR)/Sqbase.h $(CCHDR)/Sqmat.h $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Trimat.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# LU-decomposition
libinalg.a(Lu.o): $(CCHDR)/Lu.h $(CCSRC)/Lu.c++ $(CCHDR)/Sqmat.h $(CCHDR)/Vector.h $(CCTMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Lu.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# QL-decomposition and symmetric diagonalisation
libinalg.a(Ql.o): $(CCSRC)/Ql.c++ $(CCHDR)/Ql.h $(CCHDR)/Trimat.h $(CCHDR)/Sqmat.h $(CCHDR)/Vector.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Ql.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Safe division and hypot()
libinalg.a(Safety.o): $(CCSRC)/Safety.c++ $(CCHDR)/Safety.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Safety.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%
    	
# Real symmetric matrix diagonalisation: all pos. eigenvectors
libinalg.a(Rsmdiag.o): $(CCSRC)/Rsmdiag.c++ $(CCHDR)/Rsmdiag.h \
		$(CCHDR)/Trimat.h $(CCHDR)/Sqmat.h \
		$(CCHDR)/Vector.h $(CCHDR)/Safety.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Rsmdiag.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Singular value decomposition
libinalg.a(Svd.o): $(CCSRC)/Svd.c++ $(CCHDR)/Svd.h $(CCHDR)/Vmutils.h $(CCHDR)/Trimat.h \
		$(CCHDR)/Sqmat.h $(CCHDR)/Matrix.h $(CCHDR)/Vector.h $(CCHDR)/Safety.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Svd.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Vectors
libinalg.a(Vector.o): $(CCSRC)/Vector.c++ $(CCHDR)/Vector.h 
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Vector.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Other vector/matrix utilities
libinalg.a(Vmutils.o): $(CCSRC)/Vmutils.c++ $(CCHDR)/Vmutils.h $(CCHDR)/Trimat.h $(CCHDR)/Sqmat.h \
		$(CCHDR)/Matrix.h $(CCHDR)/Matbase.h $(CCHDR)/Vector.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Vmutils.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# archive
LIBINALG_SRC = $(CCSRC)/Matbase.c++ $(CCHDR)/Matbase.h \
		$(CCSRC)/Rectbase.c++ $(CCHDR)/Rectbase.h \
		$(CCSRC)/Sqbase.c++ $(CCHDR)/Sqbase.h \
		$(CCSRC)/Matrix.c++ $(CCHDR)/Matrix.h \
		$(CCSRC)/Sqmat.c++ $(CCHDR)/Sqmat.h \
		$(CCSRC)/Trimat.c++ $(CCHDR)/Trimat.h \
		$(CCSRC)/Lu.c++ $(CCHDR)/Lu.h \
		$(CCSRC)/Ql.c++ $(CCHDR)/Ql.h \
		$(CCSRC)/Rsmdiag.c++ $(CCHDR)/Rsmdiag.h \
		$(CCSRC)/Safety.c++ $(CCHDR)/Safety.h \
		$(CCSRC)/Svd.c++ $(CCHDR)/Svd.h \
		$(CCSRC)/Vector.c++ $(CCHDR)/Vector.h \
		$(CCSRC)/Vmutils.c++ $(CCHDR)/Vmutils.h

libinalg.tar: $(LIBINALG_SRC)
	tar cvf $@ $(LIBINALG_SRC)

# ---- libccstat.a ----

# The "CC statistics" library. Depends on -linalg
libccstat.a: libccstat.a(Cdf.o) libccstat.a(Pmest.o) \
		libccstat.a(Stat2.o)
	rm $?

# DSO version
# libccstat.so: libccstat.a libinalg.a
#	cd $(ABI); $(CXX) $(SOFLAGS) -shared -rdata_shared -no_unresolved -all $(@F:.so=.a) -o $(@F) -none libinalg.a

# Cumulative distribution function approximation
libccstat.a(Cdf.o): $(CCSRC)/Cdf.c++ $(CCHDR)/Cdf.h $(CCHDR)/Vector.h $(CCTMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Cdf.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# Parameter estimation
libccstat.a(Pmest.o): $(CCSRC)/Pmest.c++ $(CCHDR)/Pmest.h \
		$(CCHDR)/Vector.h $(CCHDR)/Matrix.h $(CCHDR)/Trimat.h $(CCHDR)/Vmutils.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Pmest.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Two-way statistics
libccstat.a(Stat2.o): $(CCSRC)/Stat2.c++ $(CCHDR)/Stat2.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Stat2.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# archive
LIBCCSTAT_SRC = $(CCSRC)/Cdf.c++ $(CCHDR)/Cdf.h \
		$(CCSRC)/Pmest.c++ $(CCHDR)/Pmest.h \
		$(CCSRC)/Stat2.c++ $(CCHDR)/Stat2.h

libccstat.tar: $(LIBCCSTAT_SRC)
	tar cvf $@ $(LIBCCSTAT_SRC)

# ---- libccutils.a ----

# The "ccutils" library (for miscellaneous utilities).
libccutils.a: libccutils.a(Bits.o) libccutils.a(Cplx.o) \
		libccutils.a(Spl.o) libccutils.a(String.o)
	rm $?

# DSO version
# libccutils.so: libccutils.a
#	cd $(ABI); $(CXX) $(SOFLAGS) -shared -rdata_shared -all $(@F:.so=.a) -o $(@F)

# Bit vectors
libccutils.a(Bits.o): $(CCSRC)/Bits.c++ $(CCHDR)/Bits.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Bits.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Complex arithmetics
libccutils.a(Cplx.o): $(CCSRC)/Cplx.c++ $(CCHDR)/Cplx.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/Cplx.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# Splines
libccutils.a(Spl.o): $(CCSRC)/Spl.c++ $(CCHDR)/Spl.h $(CCTMPLHDR)/Array.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Spl.c++ -o $%
	$(TMPLINK) $%
	$(AR) $(ARFLAGS) $@ $%

# strings
libccutils.a(String.o): $(CCSRC)/String.c++ $(CCHDR)/String.h
	$(CXX) $(CCFLAGS) -c $(CCSRC)/String.c++ -o $%
	$(AR) $(ARFLAGS) $@ $%

# archive
LIBCCUTILS_SRC = $(CCSRC)/Bits.c++ $(CCHDR)/Bits.h \
		$(CCSRC)/Cplx.c++ $(CCHDR)/Cplx.h \
		$(CCSRC)/Spl.c++ $(CCHDR)/Spl.h \
		$(CCSRC)/String.c++ $(CCHDR)/String.h

libccutils.tar: $(LIBCCUTILS_SRC)
	tar cvf $@ $(LIBCCUTILS_SRC)

#
# ---- C++ NON-LIBRARY MODULES ----
#

# Best RMSD rotation in high-dimensional spaces
Hirot.o: $(CCSRC)/Hirot.c++ $(CCHDR)/Hirot.h $(CCHDR)/Vector.h $(CCHDR)/Sqmat.h \
		$(CCHDR)/Lu.h $(CCHDR)/Svd.h $(CCHDR)/Points.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Hirot.c++ -o $@

# Maskable array of vectors
Points.o: $(CCSRC)/Points.c++ $(CCHDR)/Points.h $(CCTMPLHDR)/Array.h $(CCHDR)/Vector.h \
		$(CCHDR)/Bits.h $(CCHDR)/Matrix.h $(CCHDR)/Sqmat.h $(CCHDR)/Trimat.h
	$(CXX) $(CCFLAGS) $(TMPLOPTS) -c $(CCSRC)/Points.c++ -o $@

# archive
CCMODS_SRC = $(CCSRC)/Hirot.c++ $(CCHDR)/Hirot.h \
		$(CCSRC)/Points.c++ $(CCHDR)/Points.h 

ccmods.tar: $(CCMODS_SRC)
	tar cvf $@ $(CCMODS_SRC)


# ---- Template archives ----

TMPL_DECLS = $(CCTMPLHDR)/Array.h $(CCTMPLHDR)/List1.h $(CCTMPLHDR)/Maskarr.h \
		$(CCTMPLHDR)/Smartptr.h
TMPL_DEFS = $(TMPL_DECLS:.h=.c++)

templates.tar: $(TMPL_DECLS) $(TMPL_DEFS)
	tar cvf $@ $(TMPL_DECLS) $(TMPL_DEFS)

# =====================================================
