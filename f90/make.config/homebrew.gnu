# -*- Makefile -*-

#
# Setup file for GNU compiler installed with homebrew (http://brew.sh) on macOS
#
# This file is part of the JAMS Makefile system, distributed under the MIT License.
#
# Copyright (c) 2011-2019 Matthias Cuntz - mc (at) macu (dot) de
#

# Paths
GNUDIR := /usr/local
GNULIB := $(GNUDIR)/lib
GNUBIN := $(GNUDIR)/bin

# Compiling
F90 := $(GNUBIN)/gfortran
FC  := $(F90)
CC  := /usr/bin/cc  # clang, same as /usr/bin/gcc
CPP := /usr/bin/cpp # could be 'gcc -E -cpp' on Linux but does not work on Mac
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
FCFLAGS  += -ffixed-form -ffixed-line-length-132
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

# IMSL
IMSLDIR  :=
IMSLFLAG := #-limsl -limslscalar -limsllapack -limslblas -limsls_err -limslmpistub -limslsuperlu
IMSLDEF  := #-DIMSL

# MKL
INTEL  :=
MKLDIR := # $(INTEL)/mkl
ifeq ($(openmp),true)
    MKLFLAG := -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread #-lpthread
else
    MKLFLAG := -lmkl_intel_lp64 -lmkl_core -lmkl_sequential #-lpthread
endif
MKLDEF := -DMKL
MKL95DIR  :=
MKL95LIB  := # $(MKL95DIR)/lib
MKL95INC  := # $(MKL95DIR)/include/intel64/lp64
MKL95FLAG := -lmkl_lapack95_lp64
MKL95DEF  := -DMKL95

# NETCDF
ifeq ($(netcdf),netcdf3)
    NCDIR  :=
    NCFLAG := -lnetcdff -lnetcdf
    NCDEF  := -DNETCDF -DNETCDF3
else
    NCDIR    := /usr/local
    NCFLAG   := -lnetcdf
    NCDEF    := -DNETCDF
    NCFDIR   :=
    NCFFLAG  := -lnetcdff
    HDF5LIB  := /usr/local/lib
    HDF5FLAG := -lhdf5_hl -lhdf5
    SZLIB    := /usr/local/lib
    SZFLAG   := -lsz
    CURLLIB  := /usr/lib
    CURLFLAG := -lcurl
    ZFLAG    := -lz
endif

# PROJ
PROJ4DIR  := /usr/local/
PROJ4FLAG := -lproj
FPROJDIR  :=
FPROJLIB  :=
FPROJFLAG := #-lfproj4 $(FPROJLIB)/proj4.o
FPROJDEF  := #-DFPROJ

# LAPACK
LAPACKDIR  :=
LAPACKFLAG := -framework Accelerate
LAPACKDEF  := -DLAPACK

# MPI
OPENMPIDIR  :=
OPENMPIDEF  := -DMPI

# Documentation
DOXYGENDIR :=
DOTDIR     := /usr/local/bin
TEXDIR     := /Library/TeX/texbin
PERLDIR    := /usr/bin
