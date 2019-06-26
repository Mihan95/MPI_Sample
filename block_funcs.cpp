#include "prot_and_classes.h"
#include <iostream>
#include <math.h>
#include "funcs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_OUTPUT_SIZE 7

using namespace std;

double *filling(string_ind *global_index, int size, int rank, int m, int s, int l, double *A, double (*p)(int i, int j, int i1, int j1, int m),
                int process_last_row)
{
    int i  = 0;
    int j  = 0;
    int i1 = 0;
    int j1 = 0;

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    for (i = rank; i < l; i += size)
    {
        for (j = 0; j < l; j++)
        {
            for (i1 = 0; i1 < m; i1++)
                for (j1 = 0; j1 < m; j1++)
                    A[global_index[i].local_index * row + j * col_m + (i1 * m + j1)]
                            = (*p)(i, j, i1, j1, m);
        }


            for (i1 = 0; i1 < m; i1++)
                for (j1 = 0; j1 < s; j1++)
                    A[global_index[i].local_index * row + l * col_m + (i1 * s + j1)]
                            = (*p)(i, j, i1, j1, m);
    }

    if (rank == process_last_row)
    {
        for (j = 0; j < l; j++)
        {
            for (i1 = 0; i1 < s; i1++)
                for (j1 = 0; j1 < m; j1++)
                    A[global_index[l].local_index * row + j * col_sm + (i1 * m + j1)]
                            = (*p)(l, j, i1, j1, m);
        }

            for ( i1 = 0; i1 < s; i1++)
                for ( j1 = 0; j1 < s; j1++)
                    A[global_index[l].local_index * row + l * col_sm + (i1 * s + j1)] = (*p)(l, l, i1, j1, m);
    }
    return A;
}

int ReadMatrix(double *A, int m, int l, int s, const char *name, int rank, int size,
               int process_last_row, int nrow)
{
    FILE *fin = 0;

    double *B = A;

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    int dest = 0;

    int err = 0, global = 0;
    int buf_len = row * (l / size) + row;

    MPI_Status status;

    if (rank == process_last_row)
    {
        fin = fopen(name, "r");
        if (fin == 0) err = 1;
    }

    MPI_Allreduce(&err, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (global) return -1;
    if (rank == process_last_row)
    {
        memset(A, 0, buf_len * sizeof(double));

        for (int i = 0; i < l; i++)
        {
            for (int k = 0; k < m; k++)
            {
                for (int j = 0; j < l; j++)
                {
                    for (int r = 0; r < m; r++)
                    {
                          if (file_read(B[0 + j * col_m + (k * m + r)], fin) == -1)
                            return -1;
                    }
                }

                ///Остаточные столбцы
                for (int r = 0; r < s; r++)
                {
                    if (file_read(B[0 + l * col_m + (k * s + r)], fin) == -1)
                        return -1;
                }

            }

            dest = (i % size);
            if (dest != process_last_row)
                MPI_Send(B, row, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            else B += row;
        }

        ///Остаточные строки
        for (int k = 0; k < s; k++)
        {
            for (int j = 0; j < l; j++)
            {
                for (int r = 0; r < m; r++)
                {
                    if (file_read(B[0 + j * col_sm + (k * m + r)], fin) == -1)
                        return -1;
                }
            }

            ///Остаточные столбцы
            for (int r = 0; r < s; r++)
            {
                if (file_read(B[0 + l * col_sm + (k * s + r)], fin) == -1)
                    return -1;
            }
        }
    }
    else if (rank < size)
    {
        for (int i = 0; i < nrow; i++)
        {
            MPI_Recv(B, row, MPI_DOUBLE, process_last_row, 0, MPI_COMM_WORLD, &status);
            B += row;
        }
    }

    err = 0; global = 0;
    if (rank == process_last_row)
    {
        if (ferror(fin)) err = 1;
        fclose(fin);
    }

    MPI_Allreduce(&err, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (global) return -2;
    return 0;
}

int block_filling(int argc, char *argv[], double *A, int m, int l, int s, int &inputMode, int rank,
                  int size, string_ind *global_index, int process_last_row, int nrow)  ///Заполнение матрицы
{
    if (argc == 4)
    {
        return ReadMatrix(A, m, l, s, argv[3], rank, size, process_last_row, nrow);
    }
    if (argc == 3)
    {
      if (inputMode == -1)
      {
        /*cout << "Select way to assignment the matrix: \n"
        "1. Gilbert matrix\n"
        "2. Single matrix\n"
        "3. Matrix of difference\n";
        cin >> inputMode;
        cout << endl;*/
        inputMode = 3;
      }
        switch (inputMode)
        {
        case (1) :
            A = filling(global_index, size, rank, m, s, l, A, gilbert, process_last_row);
            break;
        case(2) :
            A = filling(global_index, size, rank, m, s, l, A, single, process_last_row);
            break;
        case (3) :
            A = filling(global_index, size, rank, m, s, l, A, difference, process_last_row);
            break;
        }
    }
    return 1;
}

int GetProcess(int global_i, int size)
{
    return (global_i % size);
}

void BlockPrintInFile(const double *A, int rank, int m, int l, int s, int size,
                      const string_ind *global_index)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    const char *name = "matrix.txt";
    FILE *fout;

    for (int i = 0; i < l; i++)
    {
        if (GetProcess(i, size) == rank)
        {
            fout = fopen(name, "a");
            if (fout)
            {
                for (int k = 0; k < m; k++)
                {
                    for (int j = 0; j < l; j++)
                        for (int r = 0; r < m; r++)
                            fprintf(fout, "%7.8f ", A[global_index[i].local_index * row + j * col_m + (k * m + r)]);

                        for (int r = 0; r < s; r++)
                            fprintf(fout, "%7.8f ", A[global_index[i].local_index * row + l * col_m + (k * s + r)]);

                     fprintf(fout, "\n");

                }
            }
            fclose(fout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (GetProcess(l, size) == rank)
    {
        fout = fopen(name, "a");
        if (fout)
        {
            for (int k = 0; k < s; k++)
            {
                for (int j = 0; j < l; j++)
                    for (int r = 0; r < m; r++)
                        fprintf(fout, "%7.8f ", A[global_index[l].local_index  * row + j * col_sm + (k * m + r)]);

                    for (int r = 0; r < s; r++)
                        fprintf(fout, "%7.8f ", A[global_index[l].local_index * row + l * col_sm + (k * s + r)]);

                fprintf(fout, "\n");
            }
        }
        fclose(fout);
     }
        MPI_Barrier(MPI_COMM_WORLD);
}

void BlockPrintInMonitor(double *A, int rank, int m, int l, int s, const string_ind *global_index, double *PrintBuf, int n)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    if (n > MAX_OUTPUT_SIZE)
    {
        l = MAX_OUTPUT_SIZE / m;
        s = MAX_OUTPUT_SIZE % m;
    }

    int p;
    MPI_Status status;
    for (int i = 0; i < l; i++)
    {
        p = global_index[i].process;
        if ((rank == 0) && (p == 0))
        {
            for (int k = 0; k < m; k++)
            {
                for (int j = 0; j < l; j++)
                    for (int r = 0; r < m; r++)
                        printf("%7.8f ", A[global_index[i].local_index * row + j * col_m + (k * m + r)]);

                    for (int r = 0; r < s; r++)
                        printf("%7.8f ", A[global_index[i].local_index * row + l * col_m + (k * s + r)]);

                 printf("\n");
            }
        }
        else if (rank == 0)
             {
                 MPI_Recv(PrintBuf, row, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
                 for (int k = 0; k < m; k++)
                 {
                     for (int j = 0; j < l; j++)
                         for (int r = 0; r < m; r++)
                             printf("%7.8f ", PrintBuf[0 + j * col_m + (k * m + r)]);

                         for (int r = 0; r < s; r++)
                             printf("%7.8f ", PrintBuf[0 + l * col_m + (k * s + r)]);

                      printf("\n");
                 }
             }
            else if (rank == p)
                 {
                    MPI_Send(A+global_index[i].local_index * row, row, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                 }
    }
    p = global_index[l].process;
    if ((rank == 0) && (p == 0))
    {
        for (int k = 0; k < s; k++)
        {
            for (int j = 0; j < l; j++)
                for (int r = 0; r < m; r++)
                    printf("%7.8f ", A[global_index[l].local_index  * row + j * col_sm + (k * m + r)]);

                for (int r = 0; r < s; r++)
                    printf("%7.8f ", A[global_index[l].local_index * row + l * col_sm + (k * s + r)]);

            printf("\n");
        }
    }
    else if (rank == 0)
         {
            MPI_Recv(PrintBuf, s * m * l + s * s, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
            for (int k = 0; k < s; k++)
            {
                for (int j = 0; j < l; j++)
                    for (int r = 0; r < m; r++)
                        printf("%7.8f ", PrintBuf[0 + j * col_sm + (k * m + r)]);

                    for (int r = 0; r < s; r++)
                        printf("%7.8f ", PrintBuf[0 + l * col_sm + (k * s + r)]);

                printf("\n");
            }
         }
        else if (rank == p)
             {
                MPI_Send(A+global_index[l].local_index*row, s * m * l + s * s, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
             }
}

void block_discrepancy(double *A, double *B, double *U, int m, int l, int s, int my_rank, int nrow,
                        int process_last_row, MPI_Comm comm_size, const string_ind *global_index_A,
                        double *bcasted_buffer, double *row_parts, int *local_index_B, int size)
{

    for (int i = 0; i < m*m*(l+1); i++)
    {
        bcasted_buffer[i] = 0.;
        row_parts[i] = 0.;
    }

    int nrow_1 = nrow - 1;
    double sum = 0.;
    double max = -1.;

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    double *bcasted_row = bcasted_buffer;

    for (int i = 0; i < m * m; i++)
    {
        U[i] = 0.;
    }

    memset(row_parts,      0, m * m * (l + 1) * sizeof(double));
    memset(bcasted_buffer, 0, m * m * (l + 1) * sizeof(double));

    for (int i = 0; i < l; i++)
    {
        if (my_rank != global_index_A[i].process) bcasted_row = bcasted_buffer;
        else bcasted_row = A + (global_index_A[i].local_index * row);
        MPI_Bcast(bcasted_row, m * m * l + m * s, MPI_DOUBLE, global_index_A[i].process, comm_size);

        for (int j = 0; j < l; j++)
        {
            for (int k = 0; k < nrow_1; k++)
            {
                multiplication(bcasted_row + (local_index_B[k] * col_m), B + (k * row + j * col_m), U, m, m, m);
                add(row_parts + j * m * m, U, m, m);
            }
            if (my_rank != process_last_row)
            {
                multiplication(bcasted_row + (local_index_B[nrow_1] * col_m), B + (nrow_1 * row + j * col_m), U, m, m, m);
                add(row_parts + j * m * m, U, m, m);
            }
            else
            {
                multiplication(bcasted_row + (local_index_B[nrow_1] * col_m), B + (nrow_1 * row + j * col_sm), U, m, s, m);
                add(row_parts + j * m * m, U, m, m);
            }

        }

        for (int k = 0; k < nrow_1; k++)
        {
            multiplication(bcasted_row + (local_index_B[k] * col_m), B + (k * row + l * col_m), U, m, m, s);
            add(row_parts + l * m * m, U, m, s);
        }
        if (my_rank != process_last_row)
        {
            multiplication(bcasted_row + (local_index_B[nrow_1] * col_m), B + (nrow_1 * row + l * col_m), U, m, m, s);
            add(row_parts + l * m * m, U, m, s);
        }
        else
        {
            multiplication(bcasted_row + (local_index_B[nrow_1] * col_m), B + (nrow_1 * row + l * col_sm), U, m, s, s);
            add(row_parts + l * m * m, U, m, s);
        }

        if (my_rank < size)
            MPI_Allreduce(row_parts, bcasted_buffer, l * m * m + s * m, MPI_DOUBLE, MPI_SUM, comm_size);

        sum = 0.;

        for (int k = 0; k < m; k++)
        {
            for (int j = 0; j < l; j++)
                for (int r = 0; r < m; r++)
                    sum += fabs (bcasted_buffer[j * col_m + (k * m + r)]);

            for (int r = 0; r < s; r++)
                sum += fabs (bcasted_buffer[l * col_m + (k * s + r)]);

            max = (sum > max) ? sum : max;
            sum = 0.;
        }

        memset(row_parts,   0, m * m * (l + 1) * sizeof(double));
    }

    if (my_rank != process_last_row) bcasted_row = bcasted_buffer;
    else bcasted_row = A + (global_index_A[l].local_index * row);
    MPI_Bcast(bcasted_row, s * m * l + s * s, MPI_DOUBLE, global_index_A[l].process, comm_size);

    for (int j = 0; j < l; j++)
    {
        for (int k = 0; k < nrow_1; k++)
        {
            multiplication(bcasted_row + (local_index_B[k] * col_sm), B + (k * row + j * col_m), U, s, m, m);
            add(row_parts + j * s * m, U, s, m);
        }
        if (my_rank != process_last_row)
        {
            multiplication(bcasted_row + (local_index_B[nrow_1] * col_sm), B + (nrow_1 * row + j * col_m), U, s, m, m);
            add(row_parts + j * s * m, U, s, m);
        }
        else
        {
            multiplication(bcasted_row + (local_index_B[nrow_1] * col_sm), B + (nrow_1 * row + j * col_sm), U, s, s, m);
            add(row_parts + j * s * m, U, s, m);
        }
    }

    for (int k = 0; k < nrow_1; k++)
    {
        multiplication(bcasted_row + (local_index_B[k] * col_sm), B + (k * row + l * col_m), U, s, m, s);
        add(row_parts + l * s * m, U, s, s);
    }
    if (my_rank != process_last_row)
    {
        multiplication(bcasted_row + (local_index_B[nrow_1] * col_sm), B + (nrow_1 * row + l * col_m), U, s, m, s);
        add(row_parts + l * s * m, U, s, s);
    }
    else
    {
        multiplication(bcasted_row + (local_index_B[nrow_1] * col_sm), B + (nrow_1 * row + l * col_sm), U, s, s, s);
        add(row_parts + l * s * m, U, s, s);
    }

    if (my_rank < size)
        MPI_Allreduce(row_parts, bcasted_buffer, l * m * s + s * s, MPI_DOUBLE, MPI_SUM, comm_size);

    sum = 0.;

    for (int k = 0; k < s; k ++)
    {
        for (int j = 0; j < l; j++)
            for (int r = 0; r < m; r++)
                sum += fabs(bcasted_buffer[j * col_sm + (k * m + r)]);

        for (int r = 0; r < s; r++)
            sum += fabs(bcasted_buffer[l * col_sm + (k * s + r)]);

        max = (sum > max) ? sum : max;
        sum = 0.;
     }


    if (my_rank == 0)
    {
        max = fabs(max - 1);
        printf("Discrepancy is %e\n\n", max);
    }
}

void pFindMax(double *A, double *U, double *E, int m, int s, int l, maximum &max, int *index_nb, int k,
              const string_ind *global_index, int my_rank)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    double sum  = 0.;
    double norm = 0.;
    bool isAllNormsZero = true;

    ///Ищем главный элемент
    for (int i = k; i < l; i++)
    {
        if (global_index[i].process != my_rank) continue;
        else
        {
            for (int j = k; j < l; j++)
            {
                for (int i1 = 0; i1 < m; i1++)
                {
                    for (int j1 = 0; j1 < m; j1++)
                    {
                        U[i1 * m + j1] = A[global_index[i].local_index * row + j * col_m + (i1 * m +j1)];
                        E[i1 * m + j1] = (i1 == j1) ? 1. : 0.;
                    }
                }

                if (Jordan_Inversion(U, E, m, index_nb) == -1)
                    continue;

                ///Считаем норму
                norm = 0.;
                for (int j2 = 0; j2 < m; j2++)
                    norm += fabs(E[j2]);

                for (int i2 = 0; i2 < m; i2++)
                {
                    sum = 0.;
                    for (int j2 = 0; j2 < m; j2++)
                        sum += fabs(E[i2 * m + j2]);

                    norm = (sum > norm) ? sum : norm;
                }

                if ((isAllNormsZero == true) && (norm > 0))
                {
                    max.min   = norm;
                    max.max_i = i;
                    max.max_j = j;
                    isAllNormsZero = false;
                }
                else if (isAllNormsZero == false)
                     {
                        if ((norm < max.min) && (norm > 0))
                        {
                            max.min   = norm;
                            max.max_i = i;
                            max.max_j = j;
                        }
                     }
            }
        }
    }
}

int blockFindMax(double *A, double *U, double *E, int m, int l, int s, int k, maximum max, int my_rank,
                 int size, int *index, int *index_nb, string_ind *global_index, int process_last_row, int nrow,
                 MPI_Op opMin, MPI_Comm comm_size)
{

    // k - шаг алгоритма
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    ///U - временная матрица, E - единичная
    int j = 0;
    max.min   = -1.;
    max.max_i = 0;
    max.max_j = 0;
    max.max = m;
    max.error = my_rank;

    if (my_rank < size)
    {
        max.isHere = 1;
        pFindMax(A, U, E, m, s, l, max, index_nb, k, global_index, my_rank);
    }
    else max.isHere = 0;

    if (max.min < 0)
        max.isHere = 0;

    maximum gen_max = {0, 0, 0, 0, 0, 1};

    if (my_rank < size){
    MPI_Allreduce(&max, &gen_max, sizeof(max), MPI_BYTE, opMin, comm_size);

    if (gen_max.min < 0)
    {
        if (my_rank == 0)
        {
            cout << "Jordan method is not possible" << endl;
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
    }

    string_ind temp;
    temp.local_index = global_index[k].local_index;
    temp.process     = global_index[k].process;

    global_index[k].local_index = global_index[gen_max.max_i].local_index;
    global_index[k].process     = global_index[gen_max.max_i].process;

    global_index[gen_max.max_i].local_index = temp.local_index;
    global_index[gen_max.max_i].process     = temp.process;

    ColumnTransposition(A, m, k, s, l, gen_max.max_j, process_last_row, my_rank, nrow);

    j = index[k];
    index[k] = index[gen_max.max_j];
    index[gen_max.max_j] = j;

    if (my_rank == global_index[k].process){
        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                U[i1 * m + j1] = A[global_index[k].local_index * row + k * col_m + (i1 * m + j1)];
                E[i1 * m + j1] = (i1 == j1) ? 1. : 0.;
            }
        }


        if (Jordan_Inversion(U, E, m, index_nb) == -1)
        {
            if (my_rank == 0)
            {
                cout << "Jordan method is not possible" << endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
    }
    }
    return 1;
}


void blockJordanInversion(double *A, double *B, double *U, double *E, double *MR, int my_rank, int size, int m,
                          int l, int s, maximum max, int *index, int *index_nb, int process_last_row,
                          string_ind *global_index, int nrow, MPI_Op opMin, string_ind *io_global_index, MPI_Comm comm_size)
{

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    double *C;
    double *D;

    int mm = m * m;
    int sm = s * m;
    int ss = s * s;

    int k_process;
    int k_local_index;
    int nrow_1 = nrow - 1;

    for (int k = 0; k < l; k++)
    {
        if (blockFindMax(A, U, E, m, l, s, k, max, my_rank, size, index, index_nb, global_index,
                     process_last_row, nrow, opMin, comm_size) == -1)
            break;

        if (max.error == -1)
            MPI_Abort(MPI_COMM_WORLD, 0);

        k_process     = global_index[k].process;
        k_local_index = global_index[k].local_index;

        if (my_rank == k_process)
        {
            C = A + k_local_index * row + k * col_m;
            D = B + k_local_index * row;
        }
        else
        {
            C = MR;
            D = C + mm * (l - k) + sm;
        }

        MPI_Bcast(C, mm * (l - k) + sm, MPI_DOUBLE, k_process, MPI_COMM_WORLD);
        MPI_Bcast(D, mm * l + sm,       MPI_DOUBLE, k_process, MPI_COMM_WORLD);

        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                E[i1 * m + j1] = (i1 == j1) ? 1 : 0;
                U[i1 * m + j1] = C[0 + (i1 * m + j1)];
            }
        }

        Jordan_Inversion(U, E, m, index_nb);

        DivideMainRow(C, E, U, l - k, 1, m, s);
        DivideMainRow(D, E, U, l,     0, m, s);

        C += mm;

        if (my_rank < size)
        {
            SubtractMainRow(A, A, C, U, k + 1, m, s, l, k, k_process, k_local_index, my_rank, nrow, process_last_row);
            SubtractMainRow(B, A, D, U, 0    , m, s, l, k, k_process, k_local_index, my_rank, nrow, process_last_row);
        }
    }

    //Работаем с последним(нижним правым) блоком

    k_process     = global_index[l].process;
    k_local_index = global_index[l].local_index;

    if (my_rank == k_process)
    {
        C = A + k_local_index * row + l * col_sm;
        D = B + k_local_index * row;
    }
    else
    {
        C = MR;
        D = C + ss;
    }

    MPI_Bcast(C, ss, MPI_DOUBLE, k_process, MPI_COMM_WORLD);
    MPI_Bcast(D, sm * l + ss, MPI_DOUBLE, k_process, MPI_COMM_WORLD);

    for (int i1 = 0; i1 < s; i1++)
    {
        for (int j1 = 0; j1 < s; j1++)
        {
            E[i1 * s + j1] = (i1 == j1) ? 1 : 0;
            U[i1 * s + j1] = C[(i1 * s + j1)];
        }
    }

    if (Jordan_Inversion(U, E, s, index_nb) == -1)
    {
        if (my_rank == 0)
        {
            cout << "Jordan method is not possible" << endl;
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
    }
    if (my_rank < size){
        for (int j = 0; j < l; j ++)
        {
           multiplication(E, D + (j * col_sm), U, s, s, m);

           for (int i1 = 0; i1 < s; i1++)
               for (int j1 = 0; j1 < m; j1++)
                   D[j * col_sm + (i1 * m + j1)] = U[i1 * m + j1];

        }

            multiplication(E, D + (l * col_sm), U, s, s, s);

            for (int i1 = 0; i1 < s; i1++)
               for (int j1 = 0; j1 < s; j1++)
                   D[l * col_sm + (i1 * s + j1)] = U[i1 * s + j1];


        for (int i = 0; i < nrow_1; i++)
        {
            for (int j = 0; j < l; j++)
            {
                multiplication(A + (i * row + l * col_m), D + (j * col_sm), U, m, s, m);
                subtract(B + (i * row + j * col_m), U, m, m);
            }
            multiplication(A + (i * row + l * col_m), D + (l * col_sm), U, m, s, s);
            subtract(B + (i * row + l * col_m), U, m, s);
        }
        if (my_rank != process_last_row)
        {
            for (int j = 0; j < l; j++)
            {
                multiplication(A + (nrow_1 * row + l * col_m), D + (j * col_sm), U, m, s, m);
                subtract(B + (nrow_1 * row + j * col_m), U, m, m);
            }
            multiplication(A + (nrow_1 * row + l * col_m), D + (l * col_sm), U, m, s, s);
            subtract(B + (nrow_1 * row + l * col_m), U, m, s);
        }
    }

    for(int i = 0; i < l; i++)
    {
        io_global_index[index[i]].local_index = global_index[i].local_index;
        io_global_index[index[i]].process = global_index[i].process;
    }
}

