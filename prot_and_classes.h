#include "mpi.h"

struct maximum
{
    int    max_i;
    int    max_j;
    double min;
    double max;
    int    error;
    int   isHere;
};

struct string_ind
{
    int local_index;
    int process;
};

class args_t
{
public:
    int thread;
    int p;
    int m;
    int l;
    int s;
    double *A;
    double *B;
    double **U;
    double **E;
    maximum *max;
    int *index;
    int *index_nb;
    double pr_time;
    args_t() : thread(0), p(0), m(0), l(0), s(0), A(NULL), B(NULL), U(NULL), E(NULL), pr_time(0) {}
};

int      block_filling        (int argc, char* argv[], double *A, int m, int l, int s, int &inputMode, int rank, int size, string_ind *global_index, int process_last_row, int nrow);
void     block_print          (double *A, int n, int m, int l, int s);
void     calc_param           (char* argv[], int &n, int &m, int &l, int &s, int &buf, int &size, int rank, bool &isWrongPrm, int &process_last_row, int &nrow);
void     block_discrepancy    (double *A, double *B, double *U, int m, int l, int s, int my_rank,
                               int nrow, int process_last_row, MPI_Comm comm_size, const string_ind *global_index_A,
                               double *bcasted_buffer, double *row_parts, int *local_index_B, int size);
void     blockJordanInversion (double *A, double *B, double *U, double *E, double *MR, int my_rank, int size,
                               int m, int l, int s, maximum max, int *index, int *index_nb, int process_last_row,
                               string_ind *global_index, int nrow, MPI_Op opMin, string_ind *io_global_index, MPI_Comm comm_size);
void     multiplication       (double *A, double *B, double *C, int m, int n, int r);
int      restore_A            (char *argv[], double **&A, int m, int l, int s, int &inputMode);
void     block_multiplication (double **A, double **B, double *U, int m, int l, int s);
double  *filling              (string_ind *global_index, int size, int rank, int m, int s, int l, double *A, double (*p)(int i, int j, int i1, int j1, int m), int process_last_row);
int      Jordan_Inversion     (double *A, double *B, int n, int *index);
void     subtract             (double *A, double *B, int n, int m);
double  *add                  (double *A, double *B, int n, int m);
int      RowPass              (double *A, int m, int l, int s, FILE *fin, int (p)(double &a, FILE *fin));
void     part_add             (double *A, double *B, int n, int s, int m);
void     RowTransposition     (double *A, int num, int m, int k, int s, int l, int t, int max_i, int p);
void     ColumnTransposition  (double *A, int m, int k, int s, int l, int max_j, int process_last_row, int my_rank, int nrow);
void     DivideMainRow        (double *A, double *E, double *U, int tail, int head, int m, int s);
void     SubtractMainRow      (double *B, double *A, double *C, double *U, int w, int m, int s, int l, int k, int k_process, int k_local_index, int my_rank, int nrow, int process_last_row);
void     CopyRow              (double *A, double *B, int m, int s, int l, int copy_index_A, int copy_index_B);
void     BlockPrintInFile     (const double *A, int rank, int m, int l, int s, int size, const string_ind *global_index);
void     BlockPrintInMonitor  (double *A, int rank, int m, int l, int s, const string_ind *global_index, double *PrintBuf, int n);
