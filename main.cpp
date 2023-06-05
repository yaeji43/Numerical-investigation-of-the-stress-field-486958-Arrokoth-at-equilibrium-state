#include "funcs.h"

int main() {

        clock_t start = clock();

        structure fun;


        fun.readInput(); //Read mesh file: 1)node, 2) element
        fun.CheckMesh();
        fun.defEquilbrium(); // compute stress distribution


        clock_t end = clock();
        cout << "running time:= " << (double) (end - start) / CLOCKS_PER_SEC << endl;

        return 0;
}
