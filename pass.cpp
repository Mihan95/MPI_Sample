#include "prot_and_classes.h"

using namespace std;

void ColumnTransposition(double *A, int m, int k, int s, int l, int max_j, int process_last_row, int my_rank, int nrow)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;
    int nrow_1 = nrow - 1;

    double tmp = 0.;
    for (int i = 0; i < nrow_1; i++)
    {
        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                tmp = A[i * row + max_j * col_m + (i1 * m +j1)];
                A[i * row + max_j * col_m + (i1 * m +j1)] = A[i * row + k * col_m +(i1 * m +j1)];
                A[i * row + k * col_m +(i1 * m +j1)] = tmp;
            }
        }
    }
    if (my_rank != process_last_row)
    {
        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                tmp = A[nrow_1 * row + max_j * col_m + (i1 * m +j1)];
                A[nrow_1 * row + max_j * col_m + (i1 * m +j1)] = A[nrow_1 * row + k * col_m +(i1 * m +j1)];
                A[nrow_1 * row + k * col_m +(i1 * m +j1)] = tmp;
            }
       }
    }
    else
    {
        for (int i1 = 0; i1 < s; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                tmp = A[nrow_1 * row + max_j * col_sm + (i1 * m + j1)];
                A[nrow_1 * row + max_j * col_sm + (i1 * m + j1)] = A[nrow_1 * row + k * col_sm + (i1 * m + j1)];
                A[nrow_1 * row + k * col_sm + (i1 * m + j1)] = tmp;
           }
       }
    }
}

void DivideMainRow(double *A, double *E, double *U, int tail, int head, int m, int s)
{
    int col_m  = m * m;

    for (int j = head; j < tail; j++)
    {
        multiplication(E, A + (0 + j * col_m), U, m, m, m);
        for (int i1 = 0; i1 < m; i1++)
            for (int j1 = 0; j1 < m; j1++)
                A[0 + j * col_m + (i1 * m + j1)] = U[i1 * m + j1];
    }

    multiplication(E, A + (0 + tail * col_m), U, m, m, s);
    for (int i1 = 0; i1 < m; i1++)
        for (int j1 = 0; j1 < s; j1++)
            A[0 + (tail) * col_m + (i1 * s + j1)] = U[i1 * s + j1];
}

void SubtractBloks(double *B, double *A, double *C, double *U, int w, int x, int m, int s, int l, int k, int i, int col)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    for (int j = w; j < l; j++)
    {
        multiplication(A + (i * row + k * col), C + ((j-w) * col_m), U, x, m, m);
        subtract(B + (i * row + j * col), U, x, m);
    }
    multiplication(A + (i * row + k * col), C + ((l-w) * col_m), U, x, m, s);
    subtract(B + (i * row + l * col), U, x, s);
}

void SubtractMainRow(double *B, double *A, double *C, double *U, int w, int m, int s, int l, int k, int k_process,
                     int k_local_index, int my_rank, int nrow, int process_last_row)
{
    int nrow_1 = nrow - 1;

    if (my_rank != k_process)
    {
        for (int i = 0; i < nrow_1; i++)
        {
            SubtractBloks(B, A, C, U, w, m, m, s, l, k, i, m * m);
        }
        if (my_rank != process_last_row)
            SubtractBloks(B, A, C, U, w, m, m, s, l, k, nrow_1, m * m);
        else SubtractBloks(B, A, C, U, w, s, m, s, l, k, nrow_1, s * m);
    }
    else
    {
        for (int i = 0; i < k_local_index; i++)
        {
            SubtractBloks(B, A, C, U, w, m, m, s, l, k, i, m * m);
        }
        if (k_local_index != nrow_1)
        {
            for (int i = k_local_index + 1; i < nrow_1; i++)
                SubtractBloks(B, A, C, U, w, m, m, s, l, k, i, m * m);
            if (my_rank != process_last_row)
                SubtractBloks(B, A, C, U, w, m, m, s, l, k, nrow_1, m * m);
            else SubtractBloks(B, A, C, U, w, s, m, s, l, k, nrow_1, s * m);
        }
        else return;
    }
}

void CopyRow(double *A, double *B, int m, int s, int l, int copy_index_A, int copy_index_B)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    for (int j = 0; j < l; j++)
        for (int i1 = 0; i1 < m; i1++)
            for (int j1 = 0; j1 < m; j1++)
                A[copy_index_A * row + j * col_m + (i1 * m + j1)] = B[copy_index_B * row + j * col_m + (i1 *m + j1)];

    for (int i1 = 0; i1 < m; i1++)
        for (int j1 = 0; j1 < s; j1++)
            A[copy_index_A * row + l * col_m + (i1 * s + j1)] = B[copy_index_B * row + l * col_m + (i1 * s + j1)];
}
