mesh: tcc_lt.geo
	gmsh tcc_lt.geo -3 -format msh2 -o tcc_lt.msh

interface: ./src/interface.f90
	gfortran ./src/interface.f90 -o ./interface
	rm *.mod

cfd: ./src/cuda_cfd.cuf
	nvfortran -O3 -fast ./src/cuda_cfd.cuf -o ./cfd
	rm *.mod

vtk: ./src/make_vtk.f90
	gfortran ./src/make_vtk.f90 -o ./vtk
	rm *.mod

run:
	rsync --verbose --progress --recursive . lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))
	ssh -i ~/.ssh/id_ed25519 lftimm@pcad.inf.ufrgs.br 'cd $(shell basename $(CURDIR)) && sbatch tupi.slurm'

fetch:
	rsync --verbose --progress --recursive --mkpath lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))/resultados ./
	rsync --verbose --progress --recursive --mkpath lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))/sw ./
	rsync --verbose --progress --recursive --mkpath lftimm@gppd-hpc.inf.ufrgs.br:$(shell basename $(CURDIR))/hw ./

profile:
	nsys profile --trace=cuda --force-overwrite=true -o kernel_trace --stats=true ./cfd

clean:
	rm -f cfd interface vtk
	rm -f *.mod
	rm -f tcc_lt*.v tcc_lt*.pr tcc_lt*.vts
	rm -f *.vtk cuda_times.csv
	rm -rf sw hw resultados
	rm -f *.out *.err
	rm -f kernel_trace.*

tar:
	tar --exclude=$(shell basename $(CURDIR)).tgz -czf $(shell basename $(CURDIR)).tgz .
