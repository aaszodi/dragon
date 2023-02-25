# Makefile for DRAGON-IV
# 27-Aug-1998. (c) Andras Aszodi

# The directory layout is as follows:-
#
# dragon4: root to the whole project <--YOU ARE HERE
# |
# +-- src++: *.c++, *.h files for the C++ modules
# |
# +-- src: *.[ch] files for the C modules
# |
# +-- test: c++, *.[ch] files for the module test programs
# |
# +-- bin: objects,libraries,executables
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
# |    |
# |    +-- startup: C shell scripts for env var setup
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

# ---- GENERAL DEFINITIONS ----

VERSION=4.18.1

# The shell
SHELL = /bin/sh

# Location of object files and executables
BIN = bin

# Location of the utilities
LIB = lib

# ---- ABI-DEPENDENT DEFINITIONS ----

# NOTE: this Makefile supports the maintenance of incompatible
# ABIs. Currently SGI o32,n32,n64,GCC/Linux,
# GCC/Solaris (SPARC and X86) are supported.
# The current workaround is to call another Makefile in 'bin' which
# carries out the compilation of individual objects and then performs
# the linking and puts the result into the appropriate subdirectory of 'bin'.
# The choice of the ABI depends on the value of the environment variable
# DRAGON_ABI, which is normally set by 'startup.csh'. This can be
# overridden manually (if compiling for all SGI ABIs, for example)
# by specifying the desired ABI on the make command line:-
#
# make "DRAGON_ABI=irix-n32" dragon
#
# will make DRAGON for the SGI N32 ABI.

# The ABI-specific macros are in the include files "./lib/Makefile.<ABI>"
# and "./bin/Makefile.<ABI>". The former contain general settings for
# all programs which use the "lib" utilities and will be included
# by the sub-Makefiles in both "lib" and "bin".  The latter include files
# contain macros for the DRAGON package only and will be included
# by the sub-Makefile in "bin".
# When porting to a new architecture, new include files like these
# should be constructed.

# ---- Explanation of macros ----

# The sub/Makefiles (one in "bin", the other in "lib") expect this
# macro to be set on the command line when invoked.
# ABI: selects the ABI directory
ABI = $(DRAGON_ABI)

# ---- PROGRAMS ----

# List the program (target) names here. 
PROGRAMS = dragon clumsy hbassign rank secmap sidech

help:
	@echo "Targets: " $(PROGRAMS)
	@echo "ABI: " $(ABI)

all: $(PROGRAMS)

$(BIN)/$(ABI)/:
	if ( [ ! -d $(BIN)/$(ABI) ] ) then mkdir -p $(BIN)/$(ABI); fi

$(PROGRAMS): 
	if ( [ ! -d $(LIB)/$(ABI) ] ) then mkdir $(LIB)/$(ABI); fi
	cd $(LIB)/$(ABI); $(MAKE) -f ../Makefile ABI='$(ABI)' all-c all-cc
	if ( [ ! -d $(BIN)/$(ABI) ] ) then mkdir $(BIN)/$(ABI); fi
	cd $(BIN)/$(ABI); $(MAKE) -f ../Makefile ABI='$(ABI)' VERSION='$(VERSION)' $@

#
# ---- MAKEFILES AND INSTALLATION SCRIPTS ----
#

# These are the various Makefiles which have something to do with DRAGON.
# Add a ./lib/Makefile.<ABI> to LIBMAKEFILES for each architecture.
MAKEFILES = ./Makefile ./bin/Makefile 
LIBMAKEFILES = ./lib/Makefile ./lib/Makefile.irix-o32 ./lib/Makefile.irix-n32 \
		./lib/Makefile.irix-n64 ./lib/Makefile.linux \
		./lib/Makefile.solaris-sparc ./lib/Makefile.solaris-i86
BINMAKEFILES = ./bin/Makefile ./bin/Makefile.irix-o32 ./bin/Makefile.irix-n32 \
		./bin/Makefile.irix-n64 ./bin/Makefile.linux \
		./bin/Makefile.solaris-sparc ./bin/Makefile.solaris-i86
INSTLOC = ./bin/startup
INSTSCRIPTS = startup-dev.csh \
	$(INSTLOC)/startup-irix-o32.csh \
	$(INSTLOC)/startup-irix-n32.csh \
	$(INSTLOC)/startup-irix-n64.csh \
	$(INSTLOC)/startup-linux.csh \
	$(INSTLOC)/startup-solaris-sparc.csh \
	$(INSTLOC)/startup-solaris-i86.csh

#
# ---- DATA FILES ----
#

# Default data file directory
DATA = ./data

# Default data files 
DATAFILES = $(DATA)/DEFAULT.aln $(DATA)/DEFAULT.acd $(DATA)/DEFAULT.pho \
$(DATA)/DEFAULT.par $(DATA)/DEFAULT.sim $(DATA)/pam120.sim $(DATA)/pam250.sim \
$(DATA)/DEFAULT.vol

# Example subdirectory
EXAMPLE = ./example

#
# ---- DOCUMENTATION ----
#

# default document directory
DOC = ./doc
MAN = $(DOC)/man1

# documentation
DOCDISTR = README.txt $(DOC)/ug.html $(DOC)/dragon.gif $(DOC)/infoflow.gif \
	$(DOC)/sgilogo.gif $(DOC)/penguin.gif $(DOC)/sun.gif \
	$(DOC)/gradproj.gif $(DOC)/movie.gif $(MAN)/dragon.1 $(MAN)/hbassign.1 $(MAN)/sidech.1
DOCFILES = $(DOCDISTR) 

#
# -------- HOUSEKEEPING ------
#

# srclib: saves the full DRAGON source except the test code in a tarfile.
# The full library source (which is not all needed) is also saved.

# srcown: only the DRAGON source is saved (no library)

# distr: makes a distribution. Saves all executables in ./bin and below,
# plus data files and selected items from the documentation.

# Tarfile names and compression program
TAR_SRCLIB  = dragon-srclib-$(VERSION).tar
TAR_SRCOWN = dragon-srcown-$(VERSION).tar
TAR_EXEC = dragon-$(VERSION)
COMPRESS = gzip
UNCOMPRESS = gunzip

# Source directories
SRC_OWN = src++ src
SRC_LIB = lib/cc/tmpl lib/cc/incl lib/cc/src lib/c/incl lib/c/src
SRCDIRS = $(SRC_OWN) $(SRC_LIB)

# Make the archives
# the GNU tar does not understand the '-' (read from stdin) argument
# the Emacs and Nedit backup files will be deleted
srclib-gnu:
	find $(SRCDIRS) \( -name '*~' -o -name '*.bck' \) -exec rm {} \;
	tar -cvBz -f $(TAR_SRCLIB).gz $(SRCDIRS) $(MAKEFILES) $(LIBMAKEFILES) $(BINMAKEFILES) $(INSTSCRIPTS) $(DATAFILES) $(EXAMPLE) $(DOCFILES)
srclib:
	find $(SRCDIRS) \( -name '*~' -o -name '*.bck' \) -exec rm {} \;
	find $(SRCDIRS) \( -name '*.[ch]' -o -name '*.c++' \) -follow \
		| tar cvBf $(TAR_SRCLIB) -
	tar rvBf $(TAR_SRCLIB) $(MAKEFILES) $(LIBMAKEFILES) $(BINMAKEFILES) $(INSTSCRIPTS) $(DATAFILES) $(EXAMPLE) $(DOCFILES)
	$(COMPRESS) $(TAR_SRCLIB)
srcown:
	find $(SRCDIRS) \( -name '*~' -o -name '*.bck' \) -exec rm {} \;
	find $(SRC_OWN) \( -name '*.[ch]' -o -name '*.c++' \) -follow \
		| tar cvBf $(TAR_SRCOWN) -
	tar rvBf $(TAR_SRCOWN) $(MAKEFILES) $(INSTSCRIPTS) $(DATAFILES) $(EXAMPLE) $(DOCFILES)
	$(COMPRESS) $(TAR_SRCOWN)

distr:
	cd $(BIN)/$(ABI); strip $(PROGRAMS)
	tar cvBf $(TAR_EXEC)-$(ABI).tar `find $(BIN)/$(ABI) -type f -perm -100`
	tar rvBf $(TAR_EXEC)-$(ABI).tar $(DATAFILES) $(EXAMPLE) $(DOCDISTR) $(INSTLOC)/startup-$(ABI).csh; $(COMPRESS) $(TAR_EXEC)-$(ABI).tar

# cleanup: removes all DRAGON object files and libraries
# corresponding to the current ABI
clean: clean-bin clean-lib 
clean-bin:
	if ( [ -d $(BIN)/$(ABI) ] ) then cd $(BIN)/$(ABI); rm -f *.[oa] ii_files/* ; fi
clean-lib:
	if ( [ -d $(LIB)/$(ABI) ] ) then cd $(LIB)/$(ABI); rm -f *.[oa] ii_files/* ; fi

