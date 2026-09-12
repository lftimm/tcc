mesh: tcc_lt.geo
	gmsh tcc_lt.geo -3 -format msh2 -o tcc_lt.msh

interface: ./src/interface.f90
	gfortran ./src/interface.f90 -o ./interface
	rm *.mod

cfd: ./src/cuda_cfd.cuf
	nvfortran ../src/cuda_cfd -o ./cfd
	rm *.mod

vtk: ./src/make_vtk.f90
	gfortran ./src/make_vtk.f90 -o ./vtks
	rm *.mod

run: ./src/submit.slurm
    rsync --verbose --progress --recursive . lftimm@gppd-hpc.inf.ufrgs.br:vento_cfd 
	ssh lftimm@gppd-hpc.inf.ufrgs.br "cd vento_cfd && sbatch ./src/submit.slurm"
