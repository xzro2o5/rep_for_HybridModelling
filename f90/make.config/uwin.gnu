# -*- Makefile -*-

#
# Setup file for GNU gfortran on Ubuntu on Windows
#
# This file is part of the JAMS Makefile system, distributed under the MIT License.
#
# Copyright (c) 2011-2019 Matthias Cuntz - mc (at) macu (dot) de
#

# Paths
GNUDIR := /usr/
GNULIB := $(GNUDIR)/lib
GNUBIN := $(GNUDIR)/bin

# Compiling
F90 := $(GNUBIN)/gfortran
FC  := $(F90)
CC  := $(GNUBIN)/gcc
CPP := /usr/bin/cpp
# GNU Fortran version >= 4.4
ifeq ($(release),debug)
    F90FLAGS += -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
    FCFLAGS  += -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
    CFLAGS   += -pedantic -Wall -W -O -g -Wno-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
else
    F90FLAGS += -O3 -Wno-aggressive-loop-optimizations
    FCFLAGS  += -O3 -Wno-aggressive-loop-optimizations
    CFLAGS   += -O3
endif
F90FLAGS += -cpp -ffree-form -ffixed-line-length-132
FCFLAGS  += -ffixed-form -ffixed-line-length-132 -x f77-cpp-input
CFLAGS   +=
MODFLAG  := -J# space significant
DEFINES  += -DGFORTRAN -DgFortran
# OpenMP
F90OMPFLAG := -fopenmp
FCOMPFLAG  := -fopenmp
COMPFLAG   := -fopenmp
LDOMPFLAG  := -fopenmp
OMPDEFINE  := -DOPENMP

# Linking
LIBS  += -L$(GNULIB)
RPATH += -Wl,-rpath,$(GNULIB)
# iLDPATH = $(GNUDIR)/lib/gcc/x86_64-unknown-linux-gnu/4.8.1:/usr/local/cloog/0.18.0-2/lib:/usr/local/isl/0.11.1-2/lib:/usr/local/mpc/1.0.1-3/lib:/usr/local/mpfr/3.1.2-2/lib:/usr/local/gmp/5.1.2-1/lib
# ifneq ($(LDPATH),)
#     LDPATH += :$(iLDPATH)
# else
#     LDPATH := $(iLDPATH)
# endif

# MKL
INTEL  :=
MKLDIR := #$(INTEL)/mkl
MKLINC := $(MKLDIR)/include/intel64/lp64
MKLLIB := $(MKLDIR)/lib/intel64
ifeq ($(openmp),true)
    MKLFLAG := -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread #-lpthread
else
    MKLFLAG := -lmkl_intel_lp64 -lmkl_core -lmkl_sequential #-lpthread
endif
MKLDEF := -DMKL
INTELLIB  := #$(INTEL)/compiler/lib/intel64
MKL95DIR  :=
MKL95FLAG := #-lmkl_blas95_lp64 -lmkl_lapack95_lp64
MKL95DEF  := #-DMKL95

# NETCDF
ifeq ($(netcdf),netcdf3)
    NCDIR  := /usr/local/netcdf-3.6.3-gfortran
    NCFLAG := -lnetcdff -lnetcdf
    NCDEF  := -DNETCDF -DNETCDF3
else
    NCDIR    := /usr/local
    NCFLAG   := -lnetcdf
    NCDEF    := -DNETCDF
    NCFDIR   := /usr/local/netcdf-fortran-4.4.4-gfortran
    NCFFLAG  := -lnetcdff
    HDF5LIB  := /usr/local/lib
    HDF5FLAG := -lhdf5_hl -lhdf5
    SZLIB    := /usr/local/lib
    SZFLAG   := -lsz
    ZLIB     := /usr/lib/x86_64-linux-gnu
    ZFLAG    := -lz
endif

# PROJ
PROJ4DIR  := /usr/local
PROJ4FLAG := -lproj
FPROJDIR  :=
FPROJFLAG := #-lfproj4 $(FPROJLIB)/proj4.o
FPROJDEF  := #-DFPROJ

# LAPACK
LAPACKDIR  := /usr
LAPACKFLAG := -lblas -llapack
LAPACKDEF  := -DLAPACK

# MPI
OPENMPIDIR := /usr/local/openmpi-2.0.3-gfortran
OPENMPIDEF := -DMPI

# Documentation
DOXYGENDIR := /usr/bin
DOTDIR     := /usr/bin
TEXDIR     := /usr/bin
PERLDIR    := /usr/bin
