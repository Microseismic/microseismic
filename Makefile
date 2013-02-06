all:
    ifort -openmp -o main  main.f90 sg_mod.f90
fast: 
	ifort -openmp -O3 -ipo -ip -march=native -o main main.f90  sg_mod.f90

debug:
	ifort -openmp -g -o main  main.f90 sg_mod.f90

clean:
	rm ./main
