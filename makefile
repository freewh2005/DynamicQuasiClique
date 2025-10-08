DynamicQC: main.cpp DynamicQuasiClique.hpp getArgs.hpp LBufferedMinHash.hpp
	g++ DynamicQuasiClique.hpp getArgs.hpp LBufferedMinHash.hpp main.cpp -O3 -o DynamicQC

clean:
	rm -f DynamicQC
