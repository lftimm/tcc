interface: 
	gfortran ../src/Interface_Gmsh_Hexa_v3.f90 -o ./interface

clean:
	rm *.b* *000.* *.con *.pro *.mod *.json inicia.par ./interface ./cfd

cuda:
	nvfortran ../src/CUDA_HEXA_IGOR_v1.cuf -o ./cfd
