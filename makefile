file_obj = lbm.o main.o
file_cpp = lbm.cpp main.cpp

run :
	make clean
	g++ -c $(file_cpp) -fopenmp
	g++ $(file_obj) -o main.exe -fopenmp -lpthread
	./main.exe > "outputLBM2.txt"

clean : 
	rm -f $(file_obj)

clean1 : 
	rm -f $(file_obj) *.vtr

allclean : 
	rm -f $(file_obj) *.vtr *.csv *.vtk *.txt

