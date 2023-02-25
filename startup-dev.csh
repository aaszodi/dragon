#!/bin/csh
#
# DRAGON-IV startup script for developers.
# TO BE DISTRIBUTED WITH THE SOURCE CODE ONLY.
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
set PVM_LOCATION = /usr/local/mm/progs/pvm/pvm3

# ---- Machine-specific settings ----

# ABI: select the application binary interface.
# Currently the following settings are supported:-
#
# Silicon Graphics
# irix-o32: "old" 32-bit SGI (IRIX 5.3 and above, MIPS-II, R4k,...,R10k)
# irix-n32: "new" 32-bit SGI (IRIX 6.2 and above, MIPS-III, R4k,...,R10k)
# irix-n64: 64-bit SGI (IRIX 6.2, 6.4 and above, MIPS-IV, R8k, R10k)
#
# IBM PC and compatibles
# linux: Linux 2.0, ELF executable format, Intel 486 and above
#
# Sun
# solaris-sparc: Solaris 2.x (SPARC processors)
# solaris-i86: Solaris 2.x (Intel 486 and above)
#
# Uncomment the line below and specify the desired architecture
# or leave it as it is and the script will figure out the best ABI.
# set ABI = 

# ==== No need to edit anything below this line. ====

# Data files for DRAGON
setenv DRAGON_DATA $DRAGON_ROOT/data

# Figure out operating system version, set optimal ABI if necessary
if ( ! $?ABI ) then
	set OSNAME = `uname -s`
	set RELEASE = `uname -r`

	switch ($OSNAME)
		case IRIX64:
			set ABI = irix-n64
			set IPNO = `uname -m`
			# use symmetric multiprocessing for
			# Origins and Octanes
			if ( $IPNO == "IP27" || $IPNO == "IP30" ) then
				setenv PVM_ARCH SGIMP64
			else
				setenv PVM_ARCH SGI64
			endif
			unset IPNO
		breaksw
		case IRIX*:
		switch ($RELEASE)
			case 6.*:
				set ABI = irix-n32
				setenv PVM_ARCH SGI6
			breaksw
			case 5.*:
				set ABI = irix-o32
				setenv PVM_ARCH SGI5
			breaksw
			case [34]\.*:
			echo "Your IRIX is release " $RELEASE ", please upgrade to 5.2"
			exit
		endsw
		breaksw
		case Linux:
			set ABI = linux
			setenv PVM_ARCH LINUX
		breaksw
		case SunOS:
			set SUN_MACH = `uname -m`
			switch ($SUN_MACH)
				case sun4*:
					set ABI = solaris-sparc
					setenv PVM_ARCH SUN4SOL2
				breaksw
				case i86pc:
					set ABI = solaris-i86
					setenv PVM_ARCH X86SOL2
				breaksw
			endsw
			unset SUN_MACH
		breaksw
		default:
		echo "Sorry, operating system" $OSNAME " is not supported"
		exit
	endsw
	unset OSNAME
	unset RELEASE
endif

# add executable subdirectory to path
set path = ( $path $DRAGON_ROOT/bin/$ABI )

# set up PVM if it is not done centrally
if ( ! $?PVM_ROOT ) then
	setenv PVM_ROOT $PVM_LOCATION
	set path = ( $path $PVM_ROOT/lib/$PVM_ARCH $PVM_ROOT/bin/$PVM_ARCH )
endif

# define ABI environment for compilation
setenv DRAGON_ABI $ABI

# unset temporaries
unset ABI
unset PVM_LOCATION
