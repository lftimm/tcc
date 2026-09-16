mesh: tcc_lt.geo
	gmsh tcc_lt.geo -3 -format msh2 -o tcc_lt.msh

interface: ./src/interface.f90
	gfortran ./src/interface.f90 -o ./interface
	rm -f *.mod

cfd: ./src/cuda_cfd.cuf
	nvfortran ./src/cuda_cfd.cuf -o ./cfd
	rm -f *.mod

vtk: ./src/make_vtk.f90
	gfortran ./src/make_vtk.f90 -o ./vtks
	rm -f *.mod

ssh:
	rsync --verbose --progress --recursive . lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))

tar:
	tar --exclude=$(shell basename $(CURDIR)).tgz -czf $(shell basename $(CURDIR)).tgz .
