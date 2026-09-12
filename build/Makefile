interface: 
	gfortran ../src/Interface_Gmsh_Hexa_v3.f90 -o ./interface
	rm *.mod

rerun:
	rm  *.con *.json 

clean: rerun
	rm *.b* *.pro *000.*

cuda:
	nvfortran ../src/CUDA_HEXA_IGOR_v1.cuf -o ./cfd
	rm *.mod

vtk:
	gfortran ../src/Post_VTK.f90 -o ./vtks

