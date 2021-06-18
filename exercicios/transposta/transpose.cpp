#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <ctime>
#include <fstream>
#include <cstring>
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
    cout << "acabou?" << endl;
}

void enviaMatriz(double *A, int w, int h, int m_pid, int m_nprocs){    
    MPI_Status  m_status;
    MPI_Request m_request;
    int to_pid = m_pid + 1;
    int from_pid = m_pid - 1;
    if (to_pid == 2)
    {
        return;
    }
    
    cout << "to_pid: " << to_pid << endl;
    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){
            int k = j * w + i;
            int a = A[k] * -1;
            MPI_Send(&A[k], sizeof(uint64_t), MPI_BYTE, to_pid, 0, MPI_COMM_WORLD);
            cout << "enviou: " << A[k] << " para m_pid: " << to_pid << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "acabou?" << endl;
}

void recebeMatriz(double *A, int w, int h, int m_pid, int m_nprocs){
    MPI_Status  m_status;
    MPI_Request m_request;
    int to_pid = m_pid + 1;
    int from_pid = m_pid - 1;
    for (int j = 0; j < h; j++){
        for (int i = 0; i < w; i++){
            int k = j * w + i;
            MPI_Recv(&A[k], sizeof(uint64_t), MPI_BYTE, from_pid, 0, MPI_COMM_WORLD, &m_status);
            cout << "recebeu: " << A[k] << " de m_pid: " << from_pid << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "acabou?" << endl;
}

void save2file_text(double *m_buffer, int w, int h, string filename){
  fstream output;
  cout << "Salvando arquivo: " << filename << endl;
  output.open(filename, fstream::out | fstream::trunc);
  for (int i = 0; i < h; i++){
    for (int j = 0; j < w; j++){
      int p = i * w + j;
      output << m_buffer[p] << " ";
    }
    output << endl;
  }

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
    m_offset = w * h;

    cout << "id: " << m_pid << " " << m_hostname << " offset " << m_pid * m_offset << endl;
    //MPI_File_open(MPI_COMM_WORLD, "/home/mzamith/Documents/Develop/C++/MPI/src/03-Matriz-load/mat.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &m_infile);
    //https://www.open-mpi.org/doc/v4.0/man3/MPI_File_open.3.php
    // printSerialized(m_buffer, w, h, m_pid, m_nprocs);
    // processaMatriz(m_buffer, w, h, m_pid, m_nprocs);
    posix_memalign(reinterpret_cast <void**>(&m_buffer)     , ALING, sizeof(double) * m_offset);
    posix_memalign(reinterpret_cast <void**>(&m_buffer_aux) , ALING, sizeof(double) * m_offset);
    // MPI_File_open(MPI_COMM_WORLD, "mat.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &m_infile);
    // MPI_File_seek(m_infile, m_pid * m_offset * sizeof(double), MPI_SEEK_SET);
    // MPI_File_read(m_infile, m_buffer, m_offset * sizeof(double), MPI_CHAR, &m_status);
    // MPI_File_close(&m_infile);

    if ((m_pid % 2) == 0) {
        for (uint64_t i = 0; i < (w * h); i++){
            m_buffer[i] = i; //static_cast<double> (( rand() % 1000));
        }
        cout << "enviar m_pid: " << m_pid << endl;
        enviaMatriz(m_buffer, w, h, m_pid, m_nprocs);
        cout << "enviei" << endl;
    } else {
        cout << "receber m_pid: " << m_pid << endl;
        recebeMatriz(m_buffer_aux, w, h, m_pid, m_nprocs);
        cout << "recebi" << endl;
        char* filename;
        string str_obj(to_string(m_pid)+"mat-out.bin");
        filename = &str_obj[0];
        save2file_text(m_buffer_aux, w, h, str_obj);
    }
    // MPI_File_open(MPI_COMM_WORLD, char_arr,   MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &m_outfile);
    // MPI_File_seek(m_outfile, m_pid * m_offset * sizeof(double), MPI_SEEK_SET);
    // MPI_File_write(m_outfile, m_buffer_aux, m_offset * sizeof(double), MPI_CHAR, MPI_STATUS_IGNORE);
    // MPI_File_close(&m_outfile);
    cout << "escrevi" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "sai?" << endl;
    
    

    // printSerialized(m_buffer, w, h, m_pid, m_nprocs);
    free(m_buffer);
    free(m_buffer_aux);
}
int  main (int ac, char **av) {
    MPI_Init(&ac,&av);
    start();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
