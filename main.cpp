#include "prot_and_classes.h"
#include <iostream>
#include <math.h>
#include "funcs.h"
#include <stdio.h>
#include "operations.h"
#include "get_time.h"

using namespace std;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    double w = 0;

    file_read(w, NULL);


    int n = 0; int m = 0;
    int l = 0; int s = 0;
    int process_last_row = 0;

    int size     = 0;
    int buf_len  = 0;
    int rank     = 0;
    int nrow     = 0;

    bool isWrongPrm = false;
    int  inputMode  = -1;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        gilbert(1, 1, 1, 1, 1);
        difference(1, 1, 1, 1, 1);
    }
    calc_param(argv, n, m, l, s, buf_len, size, rank, isWrongPrm, process_last_row, nrow);

    if ((isWrongPrm == true) && (rank == 0))
    {
        printf("Parmeters of command line are not correct, finish the program\n");
        return 0;
    }

/******************************************************************/

    double *A = new double [buf_len];
    double *B = new double [buf_len];

    string_ind *io_global_index = new string_ind [l + 1];
    string_ind *global_index    = new string_ind [l + 1];

    for (int i = 0; i < l + 1; i++){
        io_global_index[i].local_index = i / size;
        global_index[i].local_index    = i / size;
        io_global_index[i].process     = i % size;
        global_index[i].process        = i % size;
    }

    double *U        = new double [m * m];
    double *E        = new double [m * m];
    double *MR       = new double [2 * m * m * (l + 1)];
    int    *index_nb = new int    [m];
    double *PrintBuf = MR + (m * m * (l + 1));

    if (block_filling(argc, argv, A, m, l, s, inputMode, rank, size, io_global_index, process_last_row, nrow) < 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    B = filling(io_global_index, size, rank, m, s, l, B, single, process_last_row);

/******************************************************************/

    int buf = (l > size) ? l : size;
        buf = (buf > nrow) ? buf : nrow;
    int *index   = new int [buf];

    maximum max = {0, 0, 0, -1, 0, 1};

    MPI_Op opMin;
    MPI_Op_create((MPI_User_function *)myMin, 0, &opMin);
/***********************************************************************/
    MPI_Group group_world;
    MPI_Group group_size;
    MPI_Comm  comm_size;
    int *ranks = index;

    for(int i = 0; i < size; i++)
    {
        ranks[i] = i;
    }
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);
    MPI_Group_incl(group_world, size, ranks, &group_size);
    MPI_Comm_create(MPI_COMM_WORLD, group_size, &comm_size);

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < l; i++)
        index[i] = i;
/***********************************************************************/
    double t = get_full_time();
    blockJordanInversion(A, B, U, E, MR, rank, size, m, l, s, max, index, index_nb, process_last_row,
                         global_index, nrow, opMin, io_global_index, comm_size);
    t = get_full_time() - t;

    BlockPrintInMonitor(B, rank, m, l, s, io_global_index, PrintBuf, n);
    if (rank == 0)
        printf("Time is %f\n", t);

    for (int i = 0; i < l + 1; i++)
    {
        global_index[i].local_index = i / size;
        global_index[i].process     = i % size;
    }
    if (block_filling(argc, argv, A, m, l, s, inputMode, rank, size, global_index, process_last_row, nrow) < 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    MPI_Op_free(&opMin);


    int *local_index_B = index;

    for (int i = 0; i < l + 1; i++)
    {
        if (io_global_index[i].process == rank)
            local_index_B[io_global_index[i].local_index] = i;
    }

    if (rank < size)
        block_discrepancy(A, B, U, m, l, s, rank, nrow, process_last_row, comm_size, global_index, PrintBuf,
                          MR, local_index_B, size);

/******************************************************************/

    if (rank < size)
        MPI_Comm_free(&comm_size);

    delete [] U;
    delete [] E;
    delete [] index_nb;
    delete [] index;
    delete [] A;
    delete [] B;
    delete [] io_global_index;
    delete [] MR;
    delete [] global_index;

    MPI_Finalize();

    return 0;
}
