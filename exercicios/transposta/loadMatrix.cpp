#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdint>
#include <ctime>
#include <fstream>
#define ALING 64
using namespace std;
void printMatrix(uint64_t w, uint64_t h, double *A){
    cout << "----------------------------------" << endl;
    for (uint64_t j = 0; j < h; j++){
        for (uint64_t i = 0; i < w; i++){
            uint64_t k = j * w + i;
            cout << A[k] << " ";
        }//end-for (uint64_t i = 0; i < w; i++){
        cout << endl;
    }//end-for (uint64_t j = 0; j < h; j++){
}
int main (int ac, char **av){
    fstream input;
    uint64_t w = 10,
             h = 10;
    double *A = NULL;

    posix_memalign(reinterpret_cast <void**>(&A), ALING, sizeof(double) * w * h);
    input.open("mat-out.bin", fstream::binary | fstream::in );
    input.read(reinterpret_cast <char*>(A), sizeof(double) * w * h);
    input.close();

    printMatrix(w, h, A);

    free(A);
    return EXIT_SUCCESS;
}
