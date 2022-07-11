eigenPath= $(CURDIR)/
solverIntef= $(CURDIR)/Solver_Interfaces/
ts:
	g++ -std=c++17 -g -c Fwd_Euler.cpp 
	g++ -std=c++17 -o main.exe Fwd_Euler.o
	./main.exe
matrix_test_1:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o
	./main.exe

clean_MT_1:
	rm -v main_MatrixDemo_1.o
	rm -v main.exe
