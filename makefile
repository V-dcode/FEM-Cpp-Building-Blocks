eigenPath= $(CURDIR)/cd FEM
femLib= bin/
central_diff_truss:
	g++ -std=c++17 -g -c centralDifference_Truss.cpp -I$(femLib)
	g++ -std=c++17 -o main.exe centralDifference_Truss.o
	./main.exe
	python3 Time_intg_plotter_Truss.py

central_diff:
	g++ -std=c++17 -g -c centralDifferenceBar.cpp -I$(femLib)
	g++ -std=c++17 -o main.exe centralDifferenceBar.o
	# ./main.exe
	# python3 Time_intg_plotter.py
	
velo_verlet:
	g++ -std=c++17 -g -c Velo_Verlet.cpp -I$(femLib)
	g++ -std=c++17 -o main.exe Velo_Verlet.o
	./main.execd
	python3 Time_intg_plotter.py
	
ts:
	g++ -std=c++17 -g -c Fwd_Euler.cpp 
	g++ -std=c++17 -o main.exe Fwd_Euler.o
	./main.exe
	python3 Time_intg_plotter.py
	
matrix_test_1:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o
	./main.exe

clean_MT_1:
	rm -v main_MatrixDemo_1.o
	rm -v main.exe
