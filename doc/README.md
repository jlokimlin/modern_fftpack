# FFTPACK 5.1 - a FORTRAN library of fast Fourier transforms

12/06/2011

Authors:  Paul N. Swarztrauber and Richard A. Valent

## CAVEAT

------------------------------------------------------------------

FFTPACK5.1 is written in Fortran77 but it does not comply strictly 
with the Fortran standard.  We have not compiled a list of FFTPACK 
Fortran standards violations. Users whose projects require strict 
adherence to the standard should not use this package.

## Documentation

------------------------------------------------------------------

Documentation is provided in HTML format in file FFTPACK5.1.html.
Information about building the library follows below, as well as a 
synopsis of the library.


### Bugfixes from FFTPACK 5.0 to FFTPACK 5.1

------------------------------------------------------------------

The functionality of FFTPACK 5.1 is identical with the previous version 5.0.
The following bugfixes have been applied:

1) Corrected index error for high-level routines [CR]FFT[12][BFI] requiring 
   WSAVE array, and all dependent routines to which LWSAVE is passed. 
   Namely, in LWSAV = L + LOG(REAL(INT(L))) + 4 definitions, the summand term 
   LOG(REAL(INT(L))) has been corrected to LOG(REAL(INT(L)/LOG(2.)).

2) Corrected index error in routime C1FM1F at label 56, where array offset 
   was declared 1 that should have been 2. This resulted in C1FFTF transforms 
   of length N calculated incorrectly, where N is any prime .GT. 5, or N=7*M, 
   M .GE. 6.

3) Corrected RFFT2x routines by rewriting them. The backward transform 
   followed by forward now returns identity.  Input argument LENSAV must be
   be at least 
   L + 3*M + INT(LOG(REAL(L))/LOG(2.)) + 2*INT(LOG(REAL(M))/LOG(2.)) +12.
   Previously, the required value was smaller.

### Compiling the Library

------------------------------------------------------------------

It is best to make the library and test executables with a Fortran 90 
compiler.  The test programs use Fortran 90 intrnisics random_seed
and random_number.

The included Makefile is configured to build a static library on many
currently available unix and unix-like operating systems having a 
commonly available Fortran compiler. 

Examine file make.inc to see if your OS and compiler are represented.  
If they are not, you should modify file make.inc and the Makefile in 
each directory (main, src and test) according to your needs.

The source code is by default configured for single precision
real numbers.  If double precision is desired, Makefile and make.inc
must be modified with the appropriate compiler options for promoting
real to double precision as well as promoting constants to double
precision (this is often "-r8" on some, but not all, compilers).
