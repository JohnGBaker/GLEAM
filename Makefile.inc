# Modify these three as appropriate for the system

export ROOT
export LIB = ${ROOT}/lib
export INCLUDE = ${ROOT}/include
MCMC = ${ROOT}/ptmcmc
export CXX
export CFLAGS

ifeq ($(strip $(DISCOVER)), )
#laptop
  CXX = x86_64-apple-darwin13-g++-mp-4.8 -fcx-fortran-rules
  CC = x86_64-apple-darwin13-gcc-mp-4.8 -fcx-fortran-rules
  F90 = x86_64-apple-darwin13-gfortran-mp-4.8 -fcx-fortran-rules
  #Note the lens code seems faster (compiled with x86_64-apple-darwin13-gfortran-mp-4.8) with -O2 than -O3
  CFLAGS = -march=corei7 -g -lquadmath -Wuninitialized -O2  -fopenmp
  #CXXFLAGS = $(CFLAGS) -std=c++11 -g
  CXXFLAGS = $(CFLAGS) -std=c++11 -O2
  GSLDIR=/opt/local/lib
  GSLINC=/opt/local/include
  INC=${ROOT}/include
  LIB=${ROOT}/lib
else
#discover
  FFTW = /usr/local/other/fftw/3.1.2_gcc
# need to set this :
  NEWRAN = newran
  CXX =  icpc
  CC = icc
  CFLAGS = -Wall -g -I${FFTW}/include -wd981 -fPIC
endif

CFLAGS += -I${INCLUDE} -I${MCMC} 
LDFLAGS+= -L${LIB} -L${FFTW}/lib -lfftw3	

