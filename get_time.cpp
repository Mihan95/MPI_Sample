#include <sys/time.h>
#include <stdio.h>
#include <sys/resource.h>
#include "get_time.h"

double get_time()
{
    struct rusage buf;

    getrusage (RUSAGE_SELF, &buf);

    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1000000.;
}

double get_full_time()
{
    struct timeval buf;

    gettimeofday (&buf, 0);

    return buf.tv_sec + buf.tv_usec / 1000000.;
}
