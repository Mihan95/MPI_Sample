#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double gilbert(int i, int j, int i1, int j1, int m)
{
    return 1. / (1 + i1 + (i * m) + j1 + (j * m));
}

static double single(int i, int j, int i1, int j1, int m)
{
    (void) m;
    return ((i == j)&&(i1 == j1)) ? 1. : 0.;
}

static double difference(int i, int j, int i1, int j1, int m)
{
    return fabs((i1 + (i * m)) - (j1 + (j * m)));
}

static int file_read(double &a, FILE *fin)
{
    if (fin == NULL)
        return -1;
    if (!fscanf(fin, "%lf", &a))
    {
        printf("File error\n");
        fclose(fin);
        return -1;
    }
    return 1;
}













