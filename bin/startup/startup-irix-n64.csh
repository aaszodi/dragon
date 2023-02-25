#!/bin/csh
#
# DRAGON-IV startup script for SGI IRIX, n64 ABI.
# IRIX 6.2, 6.4 and above, MIPS-IV instruction set,
# R8k, R10k processors.
# WARNING: R10k O2-s do NOT support this ABI, use n32 instead.
#
# 14-Jan-1999.
#
# Source this script before using the programs.
#

# ---- Program locations ----

# DRAGON_ROOT: this is where the beast lives.
# If you want to make the programs available system-wide,
# a good setting would be /usr/local/dragon4.
setenv DRAGON_ROOT $HOME/dragon4

# PVM_LOCATION: you need this variable to tell DRAGON
# where it can find PVM, the Parallel Virtual Machine.
# If the PVM setup is done centrally (ie. the env var
# PVM_ROOT is already set when this script is run)
# then you do not need to set PVM_LOCATION.
set PVM_LOCATION = /usr/local/pvm/pvm3

# ==== No need to edit anything below this line. ====

# Data files for DRAGON
setenv DRAGON_DATA $DRAGON_ROOT/data

# Application binary interface
set ABI = irix-n64

# PVM architecture
set IPNO = `uname -m`
# use symmetric multiprocessing for
# Origins and Octanes
if ( $IPNO == "IP27" || $IPNO == "IP30" ) then
	setenv PVM_ARCH SGIMP64
else
	setenv PVM_ARCH SGI64
endif
unset IPNO

# add executable subdirectory to path
set path = ( $path $DRAGON_ROOT/bin/$ABI )

# set up PVM if it is not done centrally
if ( ! $?PVM_ROOT ) then
	setenv PVM_ROOT $PVM_LOCATION
	set path = ( $path $PVM_ROOT/lib/$PVM_ARCH $PVM_ROOT/bin/$PVM_ARCH )
endif

# unset temporaries
unset ABI
unset PVM_LOCATION
