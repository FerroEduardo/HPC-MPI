#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#define ALING 64
using namespace std;
void printSerialized(double *A, int w, int h, int m_pid, int m_nprocs){
    for (int p = 0; p < m_nprocs; p++){
        if (p == m_pid){
            for (int j = 0; j < h; j++){
                for (int i = 0; i < w; i++){
                    int k = j * w + i;
                    cout << A[k] << " ";
                }//end-for (int i = 0; i < w; i++){
                cout << endl;
            }//end-for (int j = 0; j < h; j++){

        }//end-if (i == m_pid){
        MPI_Barrier(MPI_COMM_WORLD);
    }//end-for (int i = 0; i < m_nprocs; i++){
}

void processaMatriz(double *A, int w, int h, int m_pid, int m_nprocs){
    MPI_Status  m_status;
    MPI_Request m_request;
    int to_pid = m_pid + 1;
    int from_pid = m_pid - 1;
    if (m_pid % 2 == 0) {
        for (int j = 0; j < h; j++){
            for (int i = 0; i < w; i++){
                int k = j * w + i;
                int a = A[k] * -1;
                MPI_Isend(&a, sizeof(double), MPI_BYTE, to_pid, 0, MPI_COMM_WORLD, &m_request);
                cout << "enviou: " << a << endl;
                MPI_Wait(&m_request, &m_status);
            }
        }
    } else {
        for (int j = 0; j < h; j++){
            for (int i = 0; i < w; i++){
                int k = j * w + i;
                MPI_Irecv(&A[k], sizeof(double), MPI_BYTE, from_pid, 0, MPI_COMM_WORLD, &m_request);
                cout << "recebeu: " << A[k] << endl;
            }
        }
    }
    // MPI_Barrier(MPI_COMM_WORLD);
}

void enviaMatriz(double *A, int w, int h, int m_pid, int m_nprocs){    
    MPI_Status  m_status;
    MPI_Request m_request;
    int to_pid = m_pid + 1;
    int from_pid = m_pid - 1;
    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){
            int k = j * w + i;
            MPI_Send(&A[k], sizeof(double), MPI_BYTE, to_pid, 0, MPI_COMM_WORLD);
            // cout << "enviou: " << A[k] << " para m_pid: " << to_pid << endl;
        }
    }
}

void recebeMatriz(double *A, int w, int h, int m_pid, int m_nprocs){
    MPI_Status  m_status;
    MPI_Request m_request;
    int to_pid = m_pid + 1;
    int from_pid = m_pid - 1;
    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){
            int k = i * w + j;
            MPI_Recv(&A[k], sizeof(double), MPI_BYTE, from_pid, 0, MPI_COMM_WORLD, &m_status);
            // cout << "recebeu: " << A[k] << " de m_pid: " << from_pid << endl;
        }
    }
}

void lerMatrizBinaria(double *m_buffer, int w, int h, string filename) {
    fstream input;
    input.open(filename, fstream::binary | fstream::in);
    input.read(reinterpret_cast <char*>(m_buffer), sizeof(double) * w * h);
    input.close();
}

void escreverMatrizBinaria(double *m_buffer, int w, int h, string filename) {
    fstream output;
    cout << "Salvando arquivo: " << filename << endl;
    output.open(filename, fstream::binary | fstream::out | fstream::trunc);
    output.write(reinterpret_cast <const char*>(m_buffer), sizeof(double) * w * h);
    output.close();
}

void start(void){
    int m_pid,
        m_nprocs;

    int w = 10,
        h = 10,
        m_offset = 0,
        m_namelen;

    char   m_hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_File m_infile, m_outfile;
    MPI_Status m_status;
    double *m_buffer = NULL;
    double *m_buffer_aux = NULL;


    MPI_Comm_size(MPI_COMM_WORLD,&m_nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&m_pid);
    MPI_Get_processor_name(m_hostname, &m_namelen);

    if ((m_nprocs % 2) != 0) {
        cout << "não é par, encerrando..." << endl;
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
  
    // h /= m_nprocs;

    cout << "id: " << m_pid << " " << m_hostname << " offset " << m_pid * w * h << endl;
    if ((m_pid % 2) == 0) {
        posix_memalign(reinterpret_cast <void**>(&m_buffer), ALING, sizeof(double) * w * h);
        lerMatrizBinaria(m_buffer, w, h, "mat.bin");
        enviaMatriz(m_buffer, w, h, m_pid, m_nprocs);
    } else {
        posix_memalign(reinterpret_cast <void**>(&m_buffer_aux), ALING, sizeof(double) * w * h);
        recebeMatriz(m_buffer_aux, w, h, m_pid, m_nprocs);
        escreverMatrizBinaria(m_buffer_aux, w, h, "mat-out.bin");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    free(m_buffer);
    free(m_buffer_aux);
}
int  main (int ac, char **av) {
    MPI_Init(&ac,&av);
    start();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
