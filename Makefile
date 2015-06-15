ROOT=${PWD}
include Makefile.inc

default: gleam
all: gleam
quad: gleam_quad

docs:
	doxygen dox.cfg

.PHONY: clean ${LIB}/libptmcmc.a ${LIB}/libprobdist.a


gleam: gleam.cc glens.cc glens.hh cmplx_roots_sg.o mlfit.hh ${MCMC}/mcmc.hh ${LIB}/libprobdist.a  ${LIB}/libptmcmc.a
	@echo "ROOT=",${ROOT}
	${CXX} ${CFLAGS} -o gleam gleam.cc glens.cc cmplx_roots_sg.o -lgsl -L${GSLDIR} -I${GSLINC} -I${MCMC} -std=c++11 -lgfortran -lprobdist -lptmcmc -L${LIB} 

gleam_quad: gleam.cc glens.cc glens.hh cmplx_roots_sg_quad.o mlfit.hh ${MCMC}/mcmc.hh ${LIB}/libprobdist.a  ${LIB}/libptmcmc.a
	${CXX} ${CFLAGS} -o gleam_quad gleam.cc glens.cc cmplx_roots_sg_quad.o -lgsl -L${GSLDIR} -I${GSLINC} -I${MCMC} -std=c++11 -lgfortran -lprobdist -lptmcmc -L${LIB} -DUSE_KIND_16 

cmplx_roots_sg.o: cmplx_roots_sg.f90
	${F90} ${CFLAGS} -c cmplx_roots_sg.f90

cmplx_roots_sg_quad.o: cmplx_roots_sg_quad.f90
	${F90} ${CFLAGS} -c cmplx_roots_sg_quad.f90

${LIB}/libptmcmc.a: ${LIB} 
	@echo "Descending to ptMCMC"
	${MAKE} -C ${MCMC}  ${MFLAGS} ${LIB}/libptmcmc.a

${LIB}/libprobdist.a: ${LIB} ${INCLUDE}
	@echo "Descending to ptMCMC"
	${MAKE} -C ${MCMC}  ${MFLAGS} ${LIB}/libprobdist.a

${LIB}:
	mkdir ${LIB}

${INCLUDE}:
	mkdir ${INCLUDE}

clean:
	rm -f *.o gleam gleam_quad
	rm lib/*.a
	${MAKE} -C ptMCMC clean
