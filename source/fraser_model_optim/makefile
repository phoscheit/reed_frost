household: fraser_model.cpp
	clang++ -std=c++11 -I /usr/local/include/eigen3 -I /usr/local/include/boost \
	-lnlopt -lm -O3 -march=native -ffast-math fraser_model.cpp -o household

simulator: simulator.cpp
	clang++ -std=c++11 -I /usr/local/include/eigen3 -I /usr/local/include/boost \
	-lnlopt -lm -O3 -march=native -ffast-math simulator.cpp -o simulator