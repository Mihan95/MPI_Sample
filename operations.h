void myMin(maximum *a, maximum *b, int *len, MPI_Datatype *datatype)
{
  (void) *len;
    for(int i = 0; i < 1; i++)
    {
        if ((a[i].min < 0) && (b[i].min < 0))
        {
            continue;
        }
        else if ((a[i].min < 0) && (b[i].min > 0))
                continue;
             else if ((a[i].min > 0) && (b[i].min < 0))
                  {
                        b[i].min   = a[i].min;
                        b[i].max_i = a[i].max_i;
                        b[i].max_j = a[i].max_j;
                        continue;
                  }

        if (a[i].min < b[i].min)
        {
            if (a[i].isHere == 1)
            {
                b[i].min   = a[i].min;
                b[i].max_i = a[i].max_i;
                b[i].max_j = a[i].max_j;
            }
            else continue;

        }
        else if (fabs(a[i].min - b[i].min) < 1e-15)
             {
                if (a[i].max_i * b[i].max + a[i].max_j > b[i].max_i * b[i].max + b[i].max_j)
                {
                    if (a[i].isHere == 1)
                    {
                        b[i].max_i = a[i].max_i;
                        b[i].max_j = a[i].max_j;
                    }
                    else continue;
                }
             }
    }

    (void) datatype;
}
