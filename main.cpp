#include <iostream>
#include <random>
#include "DynamicQuasiClique.hpp"
using namespace std;

int main(int argc, char * argv[]) {
    argsController * aC=new argsController(argc, argv);
    read_graph(*aC);
    get_res();
    return 0;
}