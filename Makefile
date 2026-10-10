mesh: tcc_lt.geo
	gmsh tcc_lt.geo -3 -format msh2 -o tcc_lt.msh

interface: ./src/interface.f90
	gfortran ./src/interface.f90 -o ./interface
	rm *.mod

cfd: ./src/cuda_cfd.cuf
	nvfortran -O3 ./src/cuda_cfd.cuf -o ./cfd
	rm *.mod

run_sanitizer: cfd
	compute-sanitizer --tool memcheck --report-api-errors explicit --leak-check full --language fortran --log-file ./sw/compute_sanitizer.txt ./cfd

vtk: ./src/make_vtk.f90
	gfortran ./src/make_vtk.f90 -o ./vtk
	rm *.mod

ssh:
	rsync --verbose --progress --recursive . lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))
	ssh -i ~/.ssh/id_ed25519 lftimm@pcad.inf.ufrgs.br 'cd $(shell basename $(CURDIR)) && sbatch tupi.slurm'
	mkdir -p resultados
	rsync --verbose --progress --recursive "lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))/tcc_lt[0-9][0-9][1-9]*" ./resultados

tar:
	tar --exclude=$(shell basename $(CURDIR)).tgz -czf $(shell basename $(CURDIR)).tgz .
