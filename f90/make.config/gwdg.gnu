# -*- Makefile -*-

#
# Setup file for GNU gfortran 6.3 on explor.univ-lorraine.fr
#
# This file is part of the JAMS Makefile system, distributed under the MIT License.
#
# Copyright (c) 2011-2019 Matthias Cuntz - mc (at) macu (dot) de
#

# Paths
GNUDIR := /cm/shared/modulefiles/gcc/9.2.0
GNULIB := $(GNUDIR)/install/lib64
GNUBIN := /usr/product/applsw/gcc/9.2.0/install/bin
#$(GNUDIR)/install/bin

# Compiling
F90 := $(GNUBIN)/gfortran
FC  := $(F90)
CC  := $(GNUBIN)/gcc
CXX := $(GNUBIN)/g++
CPP := /usr/bin/cpp
# GNU Fortran version >= 4.4
ifeq ($(release),debug)
    F90FLAGS += -pedantic-errors -Wall -static -W -O -g -Wno-maybe-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
    FCFLAGS  += -pedantic-errors -Wall -static -W -O -g -Wno-maybe-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
    CFLAGS   += -pedantic -Wall -W -O -g -Wno-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
    CXXFLAGS += -pedantic -Wall -W -O -g -Wno-uninitialized # -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal
else
    F90FLAGS += -O3 -Wno-aggressive-loop-optimizations
    FCFLAGS  += -O3 -Wno-aggressive-loop-optimizations
    CFLAGS   += -O3
    CXXFLAGS += -O3
endif
F90FLAGS += -cpp -ffree-form -ffixed-line-length-132
FCFLAGS  += -ffixed-form -ffixed-line-length-132 -x f77-cpp-input
CFLAGS   +=
CXXFLAGS +=
MODFLAG  := -J# space significant
DEFINES  += -D__GFORTRAN__ -D__gFortran__
# OpenMP
F90OMPFLAG := -fopenmp
FCOMPFLAG  := -fopenmp
COMPFLAG   := -fopenmp
CXXOMPFLAG := -fopenmp
LDOMPFLAG  := -fopenmp
OMPDEFINE  := -D__OPENMP__

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
INTEL  := /opt/soft/all/intel/compilers_and_libraries_2018.0.128/linux
MKLDIR := $(INTEL)/mkl
MKLINC := $(MKLDIR)/include
MKLLIB := $(MKLDIR)/lib/intel64_lin
ifeq ($(openmp),true)
    MKLFLAG := -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread #-lpthread
else
    MKLFLAG := -lmkl_intel_lp64 -lmkl_core -lmkl_sequential #-lpthread
endif
MKLDEF := -D__MKL__
INTELLIB  := $(INTEL)/compiler/lib/intel64_lin
MKL95DIR  :=
MKL95FLAG := #-lmkl_blas95_lp64 -lmkl_lapack95_lp64
MKL95DEF  := #-D__MKL95__

# NETCDF
ifeq ($(netcdf),netcdf3)
    NCDIR  := /home/oqx29/zzy20/local/netcdf-3.6.3-gfortran63
    NCFLAG := -lnetcdff -lnetcdf
    NCDEF  := -D__NETCDF__ -D__NETCDF3__
else
    NCDIR    := /home/oqx29/zzy20/local
    NCFLAG   := -lnetcdf
    NCDEF    := -D__NETCDF__
    NCFDIR   := /home/oqx29/zzy20/local/netcdf-fortran-4.4.4-gfortran63
    NCFFLAG  := -lnetcdff
    HDF5LIB  := /home/oqx29/zzy20/local
    HDF5FLAG := -lhdf5_hl -lhdf5
    SZLIB    := /home/oqx29/zzy20/local
    SZFLAG   := -lsz
    ZLIB     := /usr/lib64
    ZFLAG    := -lz
endif

# PROJ
PROJ4DIR  := #/usr/local/proj/4.8.0-2_gcc_4.8.1_CentOS6
PROJ4FLAG := #-lproj
FPROJDIR  :=
FPROJFLAG := #-lfproj4 $(FPROJLIB)/proj4.o
FPROJDEF  := #-D__FPROJ__

# LAPACK
LAPACKDIR  := #/usr/local/lapack/3.5.0-1_gcc_4.8.1_CentOS6
LAPACKFLAG := #-lblas -llapack
LAPACKDEF  := #-D__LAPACK__

# MPI
OPENMPIDIR := /opt/soft/hf/openmpi-2.1.1-gcc-4.9.4
OPENMPIDEF := -D__MPI__

# Documentation
DOXYGENDIR :=
DOTDIR     :=
TEXDIR     :=
PERLDIR    := /usr/bin
iiLDPATH := 
ifneq ($(LDPATH),)
    LDPATH += :$(iiLDPATH)
else
    LDPATH := $(iiLDPATH)
endif
