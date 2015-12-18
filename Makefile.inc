# Modify these three as appropriate for the system

export ROOT
export LIB = ${ROOT}/lib
export INCLUDE = ${ROOT}/include
MCMC = ${ROOT}/ptmcmc
export CXX
export CFLAGS
INC = ${ROOT}/include
LIB = ${ROOT}/lib

ifneq ($(findstring discover,$(shell hostname)),)
#discover
  #assumes: module load other/comp/gcc-5.2-sp3
  CXX = g++
  CC = gcc
  CXXFLAGS = $(CFLAGS) -std=c++11 -O3
  GSLROOT=/usr/local/other/SLES11.1/gsl/1.16/gnu-4.8.1
  CFLAGS = -O3 -lquadmath -Wuninitialized -fopenmp

else 
#laptop
  CXX = x86_64-apple-darwin13-g++-mp-4.8 -fcx-fortran-rules
  CC = x86_64-apple-darwin13-gcc-mp-4.8 -fcx-fortran-rules
  F90 = x86_64-apple-darwin13-gfortran-mp-4.8 -fcx-fortran-rules
  #Note the lens code seems faster (compiled with x86_64-apple-darwin13-gfortran-mp-4.8) with -O2 than -O3
  CFLAGS = -march=corei7 -g -lquadmath -Wuninitialized -O2  -fopenmp
  #CXXFLAGS = $(CFLAGS) -std=c++11 -g
  CXXFLAGS = $(CFLAGS) -std=c++11 -O2
  GSLROOT=/opt/local

endif

GSLDIR=$(GSLROOT)/lib
GSLINC=$(GSLROOT)/include
CFLAGS += -I${INCLUDE} -I${MCMC} 
LDFLAGS+= -L${LIB} -L${FFTW}/lib -lfftw3	

