README.txt
==========

The SGI distributions
---------------------

The DRAGON package for SGI comes in three flavors, i.e. the executables
are compiled for the three SGI application binary interfaces (ABIs).
The gzip-compressed tarfiles are called

dragon-<VERSION>-<ABI>.tar.gz

where <VERSION> is the version number like "4.18.0", and 
<ABI> is "irix-o32","irix-n32" or "irix-n64". 

Choosing the right SGI Application Binary Interface (ABI)
---------------------------------------------------------

The 'irix-o32' executables can be run on any SGI R4k/R5k/R8k/R10k system
under IRIX 5.3 and above that supports the MIPS-II instruction set.
Use this ABI if you have not upgraded from IRIX 5.3 yet for some reason.
Note that o32 will not be supported in the future.

The 'irix-n32' executables can be run on any SGI R4k/R5k/R8k/R10k system
under IRIX 6.2 that supports the MIPS-III instruction set.
I recommend that you use this ABI.

The 'irix-n64' executables can be run on SGI R8k/R10k systems under
IRIX 6.2 and above, using the MIPS-IV instruction set:
the only exceptions are O2 machines running IRIX 6.3
which do NOT support the 64-bit ABI even with an R10000 processor. The
'n64' executables were optimised for the R10000 processor
but they run on R8000-based machines as well.

There is no support for the Intel-based SGI workstations running Windows NT.

The IBM PC distribution
-----------------------

This is contained in the compressed tarfile

dragon-<VERSION>-linux.tar.gz.

DRAGON can be run on IBM PCs and compatibles 
having Intel 486, Pentium and Pentium II processors
or compatibles under Linux 2.0.x in ELF executable format.
If you have XFree86 installed and your X server
is capable of displaying 16-bit colours, then you can even enjoy
the OpenGL graphics output, thanks to Brian Paul's MESA library.
Some features do not work yet, due to bugs in GCC (version 2.8.1).

There is no DOS or Windows version of DRAGON.

The SUN distributions
---------------------

Executables are provided in the compressed tarfiles

dragon-<VERSION>-<ARCH>.tar.gz

where <ARCH> is 'solaris-sparc' for SPARC processors
or 'solaris-x86' for Intel-based systems. Unfortunately
OpenGL graphics is not supported by these executables.

Instructions for installing the DRAGON executables
--------------------------------------------------

1) Create a directory where you want DRAGON to live (say /usr/local/dragon4).

2) Copy the compressed tarfile to it. Uncompress and expand the distribution.
The executables will be expanded into the 'bin/<ARCH>' directory
where <ARCH>={irix-o32|irix-n32|irix-n64|linux|solaris-sparc|solaris-i86}.
In addition, the following subdirectories will be created:-

'data': default data files
'doc': documentation

3) Edit the 'startup-<ARCH>.csh' file in the 'bin/startup' directory
so that the variable DRAGON_ROOT is set to the
DRAGON location (/usr/local/dragon4 in our example).
You may have to edit the PVM_LOCATION variable to tell
the program where to look for PVM (if installed).
For the irix-n32 architecture, make sure that you have
the n32 ABI version of the PVM libraries.
Then source the file.
It is recommended to source 'bin/startup/startup-<ARCH>.csh' from the
.cshrc files of those users who access DRAGON regularly.

Executables
-----------

These can be found in the subdirectories of the 'bin' directory.

'dragon': the DRAGON program itself
'secmap': a program for mapping secondary structure assignments
	onto the target sequence in comparative modelling
'rank': a program that helps you to rank the models according to 
	various scores
'sidech': a partial side-chain constructing program to be used
	in comparative modelling
'hbassign': a tool for generating hydrogen bond assignments
	to be used with QUANTA
'clumsy': a simple clustering program to generate families
	from structures produced by DRAGON

PVM support
-----------

DRAGON can run in parallel using the Parallel Virtual Machine (PVM).
However, you have to obtain and install PVM yourself as it is not
part of this distribution. Try anonymous FTP from NETLIB or
a mirror site near to you.

If you want to use PVM, then set the PVM_ROOT environment variable
in 'startup-<ARCH>.csh' appropriately. Consult the PVM User Guide
for details on how to set up PVM.

Data files
----------

These live in the 'data' subdirectory. If you set up an environment
variable "DRAGON_DATA" pointing to this directory then you can refer
to these files in the DRAGON parameter file using the $DRAGON_DATA macro.
Refer to the User Guide for details.

Documentation
-------------

The complete DRAGON User Guide is the HTML document 'ug.html' in the
'doc' subdirectory. Users are strongly encouraged to read it before
attempting to use the programs. Some outdated manual pages can be found in
the 'doc/man1' subdirectory which can be moved to more convenient
locations. Manpages may be slightly out of sync with respect to the
User Guide: when in doubt, rely on the User Guide.

Bugs
----

Bug reports may be sent to me (aaszodi@nimr.mrc.ac.uk).
Reply and bug fixes are not guaranteed.

Support
-------

There is essentially none as I have discontinued active development
since I left the National Institute for Medical Research in 1996.
Sporadic updates and bug fixes may happen from time to time, though.
In general, I will not provide tutorials.

Please note that my present employer, the Novartis Forschungsinstitut GmbH
does not endorse this package in any way.

Copyright information
---------------------

Copyright (C) 1993-1997. Andras Aszodi, William R. Taylor.

The DRAGON program suite is distributed free of charge.
The copyright holders therefore undertake no warranty of any kind,
as detailed the "NO WARRANTY" section below which was taken from
the GNU General Public License:-

--------------------
  Copyright (C) 1989, 1991 Free Software Foundation, Inc.
                           675 Mass Ave, Cambridge, MA 02139, USA

                            NO WARRANTY

  BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.
--------------------
