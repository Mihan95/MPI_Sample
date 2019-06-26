#include <iostream>
#include <math.h>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS 2e-15

using namespace std;

void calc_param(char* argv[], int &n, int &m, int &l, int &s, int &buf, int &size, int rank,
                bool &isWrongParam, int &process_last_row, int &nrow)
{
    n = atoi(argv[1]);  ///Размер матрицы
    m = atoi(argv[2]);  ///Размер блоков


    if ((n == 0) || (m == 0))
    {
        isWrongParam = true;
        return;
    }

    if (n < m)
    {
        cout << "Wrong parametrs" << endl;
        exit(0);
    }
    m = (n == 2) ? 2 : m;
    ///Рассчитаем количество блоков и выделим память
    l = n / m;      ///Количество блоков размера m
    s = n % m;      ///Размер крайнего нижнего блока



    if (s == 0)
    {
        s = m;
        l--;
    }
    if (n == m)
    {
        s = 0;
        l++;
    }

    process_last_row = (l + 1) % size - 1;
    if (process_last_row < 0) process_last_row += size;

    if (l + 1  < size)
    {
        size = l + 1;
        process_last_row = l;
    }

    buf = (m * m * l + s * m) * (l / size);
    buf += m * m * l + s * m;
    /*
    if (rank < l % size)


    if (rank == l % size)
        buf += s * m * l + s * s;

    if (rank == process_last_row)
        buf = (m * m * l + s * m) * (l / size) + (m * m * l + s * m);*/

    // Количество строк в процессе
    nrow = l / size + (rank > l % size ? 0 : 1);
}

/*int find_max(double *A, double *B, int n, int k, int &max_i, int &max_j)
{
    ///A и B размера n
    ///k - шаг алгоритма

    #define A(i,j) A[(i) * n + (j)]
    #define B(i,j) B[(i) * n + (j)]
    double max = fabs(A(k,k)); max_i = k; max_j = k;

    ///Находим главный элемент
    for (int i = k; i < n; i++)
    {
        for (int j = k; j < n; j++)
        {
            if (fabs(A(i,j)) > max)
            {
                max   = fabs(A(i,j));
                max_i = i;
                max_j = j;
            }
        }
    }

    if (max < EPS)
    {
        return -1;
    }

    double temp = 0.;

    ///Переставляем строки A
    for (int j = 0; j < n; j++)
    {
        temp        = A(max_i, j);
        A(max_i, j) = A(k, j);
        A(k, j)     = temp;
    }

    ///Переставляем строки B
    for (int j = 0; j < n; j++)
    {
        temp        = B(max_i, j);
        B(max_i, j) = B(k, j);
        B(k, j)     = temp;
    }

    ///Переставляем столбцы A
    for (int i = 0; i < n; i++)
    {
        temp        = A(i, max_j);
        A(i, max_j) = A(i, k);
        A(i, k)     = temp;
    }

    #undef A
    #undef B
    return 1;
}*/

int find_max(double *A, double *B, int n, int k, int &max_i, int &max_j, int m)
{
    ///A и B размера n
    ///k - шаг алгоритма

    #define A(i,j) A[(i) * n + (j)]
    #define B(i,j) B[(i) * n + (j)]
    double max = fabs(A(k,k)); max_i = k; max_j = k;

    ///Находим главный элемент
    for (int i = k; i < n; i++)
    {
        for (int j = k; j < n; j++)
        {
            if (fabs(A(i,j)) > max)
            {
                max   = fabs(A(i,j));
                max_i = i;
                max_j = j;
            }
        }
    }

    if (max < EPS)
    {
        return -1;
    }

    double temp1 = 0., temp2 = 0., temp3 = 0., temp4 = 0.;

    ///Переставляем строки A
    for (int j = 0; j < m; j += 4)
    {
        temp1 = A(max_i, j);
        temp2 = A(max_i, j + 1);
        temp3 = A(max_i, j + 2);
        temp4 = A(max_i, j + 3);


        A(max_i, j)     = A(k, j);
        A(max_i, j + 1) = A(k, j + 1);
        A(max_i, j + 2) = A(k, j + 2);
        A(max_i, j + 3) = A(k, j + 3);

        A(k, j)     = temp1;
        A(k, j + 1) = temp2;
        A(k, j + 2) = temp3;
        A(k, j + 3) = temp4;
    }
    for (int j = m; j < n; j++)
    {
        temp1       = A(max_i, j);
        A(max_i, j) = A(k, j);
        A(k, j)     = temp1;
    }

    ///Переставляем строки B
    for (int j = 0; j < m; j += 4)
    {
        temp1 = B(max_i, j);
        temp2 = B(max_i, j + 1);
        temp3 = B(max_i, j + 2);
        temp4 = B(max_i, j + 3);

        B(max_i, j)     = B(k, j);
        B(max_i, j + 1) = B(k, j + 1);
        B(max_i, j + 2) = B(k, j + 2);
        B(max_i, j + 3) = B(k, j + 3);

        B(k, j)     = temp1;
        B(k, j + 1) = temp2;
        B(k, j + 2) = temp3;
        B(k, j + 3) = temp4;
    }
    for (int j = m; j < n; j++)
    {
        temp1       = B(max_i, j);
        B(max_i, j) = B(k, j);
        B(k, j)     = temp1;
    }

    ///Переставляем столбцы A
    for (int i = 0; i < m; i += 4)
    {
        temp1 = A(i, max_j);
        temp2 = A(i + 1, max_j);
        temp3 = A(i + 2, max_j);
        temp4 = A(i + 3, max_j);

        A(i, max_j)     = A(i, k);
        A(i + 1, max_j) = A(i + 1, k);
        A(i + 2, max_j) = A(i + 2, k);
        A(i + 3, max_j) = A(i + 3, k);

        A(i, k)     = temp1;
        A(i + 1, k) = temp2;
        A(i + 2, k) = temp3;
        A(i + 3, k) = temp4;
    }
    for (int i = m; i < n; i++)
    {
        temp1       = A(i, max_j);
        A(i, max_j) = A(i, k);
        A(i, k)     = temp1;
    }

    #undef A
    #undef B
    return 1;
}

int Jordan_Inversion(double *A, double *B, int n, int *index)
{
    #define A(i,j) A[(i) * n + (j)]
    #define B(i,j) B[(i) * n + (j)]

    int i = 0; int j = 0; int k = 0;

    int m = n - n % 4;
    int l = 0;
    int s = 0;
    int t = 0;
    int r = 0;
    int q = n - n % 3;

    for (i = 0; i < m; i += 4)
    {
        index[i]     = i;
        index[i + 1] = i + 1;
        index[i + 2] = i + 2;
        index[i + 3] = i + 3;
    }
    for (i = m; i < n; ++i)
        index[i] = i;

    double temp = 0.;
    int max_j   = 0;
    int max_i   = 0;
    int max     = 0;

    for (k = 0; k < n; k++)
    {
        max = 0;
        max = find_max(A, B, n, k, max_i, max_j, m);
        if (max == -1)
          return -1;

        j            = index[k];
        index[k]     = index[max_j];
        index[max_j] = j;

        temp = A(k,k);

        ///Делим "первую"(в данном шаге алгоритма) строку на первый элемент
        l = n - k - 1;
        s = n - l % 4;
        for (j = k + 1; j < s; j += 4)
        {
            A(k, j)     = A(k, j)     / temp;
            A(k, j + 1) = A(k, j + 1) / temp;
            A(k, j + 2) = A(k, j + 2) / temp;
            A(k, j + 3) = A(k, j + 3) / temp;
        }
        for (j = s; j < n; j++)
            A(k, j) = A(k, j) / temp;

        A(k,k) = 1;

        ///Делим "первую"(в данном шаге алгоритма) строку на первый элемент
        for (j = 0; j < m; j += 4)
        {
            B(k, j)     = B(k, j)     / temp;
            B(k, j + 1) = B(k, j + 1) / temp;
            B(k, j + 2) = B(k, j + 2) / temp;
            B(k, j + 3) = B(k, j + 3) / temp;
        }
        for (j = m; j < n; j++)
            B(k, j) = B(k, j) / temp;

        ///Вычитаем "первую" строку из нижних и верхних
        t = k - k % 3;
        r = n - l % 3;
        for (i = 0; i < t; i += 3)
        {
            for (j = k + 1; j < r; j += 3)
            {
                A(i, j)         = A(i, j)         - A(i, k)     * A(k, j);
                A(i + 1, j)     = A(i + 1, j)     - A(i + 1, k) * A(k, j);
                A(i + 2, j)     = A(i + 2, j)     - A(i + 2, k) * A(k, j);
                A(i, j + 1)     = A(i, j + 1)     - A(i, k)     * A(k, j + 1);
                A(i + 1, j + 1) = A(i + 1, j + 1) - A(i + 1, k) * A(k, j + 1);
                A(i + 2, j + 1) = A(i + 2, j + 1) - A(i + 2, k) * A(k, j + 1);
                A(i, j + 2)     = A(i, j + 2)     - A(i, k)     * A(k, j + 2);
                A(i + 1, j + 2) = A(i + 1, j + 2) - A(i + 1, k) * A(k, j + 2);
                A(i + 2, j + 2) = A(i + 2, j + 2) - A(i + 2, k) * A(k, j + 2);
            }
            for (j = r; j < n; j++)
            {
                A(i, j)     = A(i, j)     - A(i, k)     * A(k, j);
                A(i + 1, j) = A(i + 1, j) - A(i + 1, k) * A(k, j);
                A(i + 2, j) = A(i + 2, j) - A(i + 2, k) * A(k, j);
            }

        }
        for (i = t; i < k; i++)
        {
            for (j = k + 1; j < r; j += 3)
            {
                A(i, j)     = A(i, j)     - A(i, k) * A(k, j);
                A(i, j + 1) = A(i, j + 1) - A(i, k) * A(k, j + 1);
                A(i, j + 2) = A(i, j + 2) - A(i, k) * A(k, j + 2);
            }
            for (j = r; j < n; j++)
                A(i, j) = A(i, j) - A(i, k) * A(k, j);
        }

        for (i = k + 1; i < r; i += 3)
        {
            for (j = k + 1; j < r; j += 3)
            {
                A(i, j)         = A(i, j)         - A(i, k)     * A(k, j);
                A(i + 1, j)     = A(i + 1, j)     - A(i + 1, k) * A(k, j);
                A(i + 2, j)     = A(i + 2, j)     - A(i + 2, k) * A(k, j);
                A(i, j + 1)     = A(i, j + 1)     - A(i, k)     * A(k, j + 1);
                A(i + 1, j + 1) = A(i + 1, j + 1) - A(i + 1, k) * A(k, j + 1);
                A(i + 2, j + 1) = A(i + 2, j + 1) - A(i + 2, k) * A(k, j + 1);
                A(i, j + 2)     = A(i, j + 2)     - A(i, k)     * A(k, j + 2);
                A(i + 1, j + 2) = A(i + 1, j + 2) - A(i + 1, k) * A(k, j + 2);
                A(i + 2, j + 2) = A(i + 2, j + 2) - A(i + 2, k) * A(k, j + 2);
            }
            for (j = r; j < n; j++)
            {
                A(i, j)     = A(i, j)     - A(i, k)     * A(k, j);
                A(i + 1, j) = A(i + 1, j) - A(i + 1, k) * A(k, j);
                A(i + 2, j) = A(i + 2, j) - A(i + 2, k) * A(k, j);
            }

        }
        for (i = r; i < n; i++)
        {
            for (j = k + 1; j < r; j += 3)
            {
                A(i, j)     = A(i, j)     - A(i, k) * A(k, j);
                A(i, j + 1) = A(i, j + 1) - A(i, k) * A(k, j + 1);
                A(i, j + 2) = A(i, j + 2) - A(i, k) * A(k, j + 2);
            }
            for (j = r; j < n; j++)
                A(i, j) = A(i, j) - A(i, k) * A(k, j);
        }

        ///Вычитаем "первую" строку из нижних и верхних
        for (i = 0; i < t; i += 3)
        {
            for (j = 0; j < q; j += 3)
            {
                B(i, j)         = B(i, j)         - A(i, k)     * B(k, j);
                B(i + 1, j)     = B(i + 1, j)     - A(i + 1, k) * B(k, j);
                B(i + 2, j)     = B(i + 2, j)     - A(i + 2, k) * B(k, j);
                B(i, j + 1)     = B(i, j + 1)     - A(i, k)     * B(k, j + 1);
                B(i + 1, j + 1) = B(i + 1, j + 1) - A(i + 1, k) * B(k, j + 1);
                B(i + 2, j + 1) = B(i + 2, j + 1) - A(i + 2, k) * B(k, j + 1);
                B(i, j + 2)     = B(i, j + 2)     - A(i, k)     * B(k, j + 2);
                B(i + 1, j + 2) = B(i + 1, j + 2) - A(i + 1, k) * B(k, j + 2);
                B(i + 2, j + 2) = B(i + 2, j + 2) - A(i + 2, k) * B(k, j + 2);
            }
            for (j = q; j < n; j++)
            {
                B(i, j)     = B(i, j)     - A(i, k)     * B(k, j);
                B(i + 1, j) = B(i + 1, j) - A(i + 1, k) * B(k, j);
                B(i + 2, j) = B(i + 2, j) - A(i + 2, k) * B(k, j);
            }

        }
        for (i = t; i < k; i++)
        {
            for (j = 0; j < q; j += 3)
            {
                B(i, j)     = B(i, j)     - A(i, k) * B(k, j);
                B(i, j + 1) = B(i, j + 1) - A(i, k) * B(k, j + 1);
                B(i, j + 2) = B(i, j + 2) - A(i, k) * B(k, j + 2);
            }
            for (j = q; j < n; j++)
                B(i, j) = B(i, j) - A(i, k) * B(k, j);
        }

        for (i = k + 1; i < r; i += 3)
        {
            for (j = 0; j < q; j += 3)
            {
                B(i, j)         = B(i, j)         - A(i, k)     * B(k, j);
                B(i + 1, j)     = B(i + 1, j)     - A(i + 1, k) * B(k, j);
                B(i + 2, j)     = B(i + 2, j)     - A(i + 2, k) * B(k, j);
                B(i, j + 1)     = B(i, j + 1)     - A(i, k)     * B(k, j + 1);
                B(i + 1, j + 1) = B(i + 1, j + 1) - A(i + 1, k) * B(k, j + 1);
                B(i + 2, j + 1) = B(i + 2, j + 1) - A(i + 2, k) * B(k, j + 1);
                B(i, j + 2)     = B(i, j + 2)     - A(i, k)     * B(k, j + 2);
                B(i + 1, j + 2) = B(i + 1, j + 2) - A(i + 1, k) * B(k, j + 2);
                B(i + 2, j + 2) = B(i + 2, j + 2) - A(i + 2, k) * B(k, j + 2);
            }
            for (j = q; j < n; j++)
            {
                B(i, j)     = B(i, j)     - A(i, k)     * B(k, j);
                B(i + 1, j) = B(i + 1, j) - A(i + 1, k) * B(k, j);
                B(i + 2, j) = B(i + 2, j) - A(i + 2, k) * B(k, j);
            }

        }
        for (i = r; i < n; i++)
        {
            for (j = 0; j < q; j += 3)
            {
                B(i, j)     = B(i, j)     - A(i, k) * B(k, j);
                B(i, j + 1) = B(i, j + 1) - A(i, k) * B(k, j + 1);
                B(i, j + 2) = B(i, j + 2) - A(i, k) * B(k, j + 2);
            }
            for (j = q; j < n; j++)
                B(i, j) = B(i, j) - A(i, k) * B(k, j);
        }

    }

    ///Переставляем строки B обратно
    for (i = 0; i < q; i += 3)
    {
        for (j = 0; j < q; j += 3)
        {
            A(index[i],j)         = B(i,j);
            A(index[i],j + 1)     = B(i,j + 1);
            A(index[i],j + 2)     = B(i,j + 2);
            A(index[i + 1],j)     = B(i + 1,j);
            A(index[i + 1],j + 1) = B(i + 1,j + 1);
            A(index[i + 1],j + 2) = B(i + 1,j + 2);
            A(index[i + 2],j)     = B(i + 2,j);
            A(index[i + 2],j + 1) = B(i + 2,j + 1);
            A(index[i + 2],j + 2) = B(i + 2,j + 2);
        }
        for (j = q; j < n; j++)
        {
            A(index[i],j)     = B(i,j);
            A(index[i + 1],j) = B(i + 1,j);
            A(index[i + 2],j) = B(i + 2,j);
        }
    }
    for (i = q; i < n; i++)
    {
        for (j = 0; j < q; j += 3)
        {
            A(index[i],j)     = B(i,j);
            A(index[i],j + 1) = B(i,j + 1);
            A(index[i],j + 2) = B(i,j + 2);
        }
        for (j = q; j < n; j++)
        {
            A(index[i],j) = B(i,j);
        }
    }

    for (i = 0; i < q; i += 3)
    {
        for (j = 0; j < q; j += 3)
        {
            B(i,j)         = A(i,j);
            B(i,j + 1)     = A(i,j + 1);
            B(i,j + 2)     = A(i,j + 2);
            B(i + 1,j)     = A(i + 1,j);
            B(i + 1,j + 1) = A(i + 1,j + 1);
            B(i + 1,j + 2) = A(i + 1,j + 2);
            B(i + 2,j)     = A(i + 2,j);
            B(i + 2,j + 1) = A(i + 2,j + 1);
            B(i + 2,j + 2) = A(i + 2,j + 2);
        }
        for (j = q; j < n; j++)
        {
            B(i,j)     = A(i,j);
            B(i + 1,j) = A(i + 1,j);
            B(i + 2,j) = A(i + 2,j);
        }
    }
    for (i = q; i < n; i++)
    {
        for (j = 0; j < q; j += 3)
        {
            B(i,j)     = A(i,j);
            B(i,j + 1) = A(i,j + 1);
            B(i,j + 2) = A(i,j + 2);
        }
        for (j = q; j < n; j++)
        {
            B(i,j) = A(i,j);
        }
    }

    #undef A
    #undef B
    return 1;
}

/*int Jordan_Inversion(double *A, double *B, int n, int *index)
{
    #define A(i,j) A[(i) * n + (j)]
    #define B(i,j) B[(i) * n + (j)]

    int i = 0; int j = 0; int k = 0;

    for (int i = 0; i < n; ++i)
        index[i] = i;

    double temp = 0.;
    int max_j   = 0;
    int max_i   = 0;
    int max     = 0;

    for (k = 0; k < n; k++)
    {
        max = 0;
        max = find_max(A, B, n, k, max_i, max_j);
        if (max == -1)
          return -1;

        j            = index[k];
        index[k]     = index[max_j];
        index[max_j] = j;

        temp = A(k,k);

        ///Делим "первую"(в данном шаге алгоритма) строку на первый элемент
        for (j = k + 1; j < n; j++)
            A(k, j) = A(k, j) / temp;

        A(k,k) = 1;

        ///Делим "первую"(в данном шаге алгоритма) строку на первый элемент
        for (j = 0; j < n; j++)
            B(k, j) = B(k, j) / temp;

        ///Вычитаем "первую" строку из нижних и верхних
        for (i = 0; i < k; i++)
            for (j = k + 1; j < n; j++)
                A(i,j) = A(i,j) - A(i,k) * A(k,j);

        for (i = k + 1; i < n; i++)
            for (j = k + 1; j < n; j++)
                A(i,j) = A(i,j) - A(i,k) * A(k,j);

        ///Вычитаем "первую" строку из нижних и верхних
        for (i = 0; i < k; i++)
            for (j = 0; j < n; j++)
                B(i,j) = B(i,j) - A(i,k) * B(k,j);

        for (i = k+1; i < n; i++)
            for ( j = 0; j < n; j++)
                B(i,j) = B(i,j) - A(i,k)*B(k,j);

    }

    ///Переставляем строки B обратно
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            A(index[i],j) = B(i,j);

    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            B(i,j) = A(i,j);

    #undef A
    #undef B
    return 1;
}*/

void multiplication(double *A, double *B, double *C, int n, int m, int l){
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * l + (j)]
    #define C(i,j) C[(i) * l + (j)]

    int i, j, k;

    double s0, s1, s2, s3, s4, s5, s6, s7, s8;

    memset(C, 0, n * l * sizeof(double));

    int on = n - n % 3, ol = l - l % 3;

    for (i = 0; i < on; i+=3)
    {
        for (j = 0; j < ol; j+=3)
        {
        s0 = s1 = s2 = s3 = s4 = s5 = s6 = s7 = s8 = 0.;
            for (k = 0; k < m; k++)
            {
                s0 += A(i,k)   * B(k,j);
                s1 += A(i+1,k) * B(k,j);
                s2 += A(i+2,k) * B(k,j);
                s3 += A(i,k)   * B(k,j+1);
                s4 += A(i+1,k) * B(k,j+1);
                s5 += A(i+2,k) * B(k,j+1);
                s6 += A(i,k)   * B(k,j+2);
                s7 += A(i+1,k) * B(k,j+2);
                s8 += A(i+2,k) * B(k,j+2);
            }
            C(i,j)      += s0;
            C(i+1,j)    += s1;
            C(i+2,j)    += s2;
            C(i,j+1)    += s3;
            C(i+1,j+1)  += s4;
            C(i+2, j+1) += s5;
            C(i,j+2)    += s6;
            C(i+1,j+2)  += s7;
            C(i+2,j+2)  += s8;
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = ol; j < l; j++)
        {
            C(i,j) = 0;
            for (k = 0; k < m; k++)
                C(i,j) += A(i,k) * B(k,j);
        }
    }

    for (i = on; i < n; i++)
    {
        for (j = 0; j < l; j++)
        {
            C(i,j) = 0;
            for (k = 0; k < m; k++)
                C(i,j) += A(i,k) * B(k,j);
        }
    }

    #undef A
    #undef B
    #undef C
}

void subtract(double *A, double *B, int n, int m)
{
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * m + (j)]

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A(i,j) -= B(i,j);

    #undef A
    #undef B
}

void add(double *A, double *B, int n, int m)
{
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * m + (j)]

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A(i,j) += B(i,j);

    #undef A
    #undef B
}

void part_add(double *A, double *B, int n, int s, int m)
{
    m += 0;
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * s + (j)]

    for (int i = 0; i < n; i++)
        for (int j = 0; j < s; j++)
            A(i,j) += B(i,j);

    #undef A
    #undef B
}


