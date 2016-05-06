
LIB_NAME = fftpack

LIB=../lib/lib$(LIB_NAME).a

UNAMES := $(shell uname -s)

ifeq ($(UNAMES),Linux)

  PGI := $(shell pgf90 2>&1)

  ifeq ($(PGI),pgf90-Warning-No files to process)

    FC := pgf90 -module ../lib -I../lib
    CPP := pgf90 -E

  else

    FC := gfortran -Wall \
    -J../lib -I../lib
    CPP := gfortran -E
   
  endif

  MAKE := gmake
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),AIX)

  FC := xlf -qmoddir=../lib -I../lib
  CPP := xlf -d -qnoobject
  MAKE := gmake
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),SunOS)

    AR := /usr/ccs/bin/ar
    FC := /opt/SUNWspro/bin/f90 -moddir=../lib -I../lib
    CPP := /opt/SUNWspro/bin/f90 -F
    MAKE := /fs/local/bin/make

endif

ifeq ($(UNAMES),IRIX64)

    AR := /usr/bin/ar
    FC := f90 -I../lib
    CPP := f90 -E
    MAKE := /usr/local/bin/gmake

endif

ifeq ($(UNAMES),Darwin)

    AR := /usr/local/bin/ar
    FC := gfortran -O3 -Wall
    CPP := gfortran -cpp
    MAKE := /usr/local/bin/gmake

endif
