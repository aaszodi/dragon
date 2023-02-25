#!/bin/csh
#
# DRAGON-IV startup script for SGI IRIX, n32 ABI.
# IRIX 6.2 and above, MIPS-III instruction set,
# R4k, R5k, R8k, R10k processors.
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

# PVM architecture: PVM 3.10 used the 64-bit ABI
# for R8k and R10k processors and o32 for R4k processors
# and there was no separate SGIN32 architecture
# so check your setup carefully - you may need to hack around a bit.
setenv PVM_ARCH SGI6

# ==== No need to edit anything below this line. ====

# Data files for DRAGON
setenv DRAGON_DATA $DRAGON_ROOT/data

# Application binary interface
set ABI = irix-n32

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
