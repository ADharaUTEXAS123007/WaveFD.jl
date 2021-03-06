#include<stdio.h>

void
srcillum_helper_float(
    float* field_ginsu,
    size_t n,
    int    nthreads)
{
    size_t i;
#pragma omp parallel for num_threads(nthreads)
    for (i = 0; i < n; i++) {
        field_ginsu[i] = field_ginsu[i]*field_ginsu[i];
    }
}
void
srcillum_helper_double(
    double* field_ginsu,
    size_t  n,
    int     nthreads)
{
    size_t i;
#pragma omp parallel for num_threads(nthreads)
    for (i = 0; i < n; i++) {
        field_ginsu[i] = field_ginsu[i]*field_ginsu[i];
    }
}

void
illum_accumulate_float(
    float* field_ginsu_accum,
    float* field_ginsu,
    size_t n,
    int    nthreads)
{
    size_t i;
#pragma omp parallel for num_threads(nthreads)
    for (i = 0; i < n; i++) {
        field_ginsu_accum[i] += field_ginsu[i]*field_ginsu[i];
    }
}
void
illum_accumulate_double(
    double* field_ginsu_accum,
    double* field_ginsu,
    size_t n,
    int    nthreads)
{
    size_t i;
#pragma omp parallel for num_threads(nthreads)
    for (i = 0; i < n; i++) {
        field_ginsu_accum[i] += field_ginsu[i]*field_ginsu[i];
    }
}
