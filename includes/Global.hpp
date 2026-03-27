#pragma once
#include <cmath>

#if defined __MIC__
    #define VECLEN 8
    #define ALIGN  64
#elif defined __AVX512F__
    #define VECLEN 8
    #define ALIGN  64
#elif defined __AVX__
    #define VECLEN 4
    #define ALIGN  32
#else
    #define VECLEN 2
    #define ALIGN  16
#endif

#define FORCE_INLINE __attribute__((always_inline)) inline
#define ALIGNED      __attribute__((aligned(ALIGN)))

#define NDIMS 2
#define NVARS 4

template <typename T>
void setarr(int i1, int i2, T val, T* arr) {
    for (int i = i1; i < i2; i++) {
        arr[i] = val;
    }
}

FORCE_INLINE int max(int x, int y) {
    return (x >= y) ? x : y;
}

FORCE_INLINE int ceil(int x, int y) {
    return (x + y - 1) / y;
}

FORCE_INLINE double dot(const double x1[NDIMS], const double x2[NDIMS]) {
    #if   NDIMS == 2
    return x1[0]*x2[0] + x1[1]*x2[1];
    #elif NDIMS == 3
    return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
    #endif
}

FORCE_INLINE void dotv(const double x1[NDIMS][VECLEN], const double x2[NDIMS][VECLEN], double res[VECLEN]) {
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        #if   NDIMS == 2
        res[v] = x1[0][v]*x2[0][v] + x1[1][v]*x2[1][v];
        #elif NDIMS == 3
        res[v] = x1[0][v]*x2[0][v] + x1[1][v]*x2[1][v] + x1[2][v]*x2[2][v];
        #endif
    }
}

FORCE_INLINE double vabs(const double x[NDIMS]) {
    return sqrt(dot(x,x));
}

FORCE_INLINE void vabsv(const double x[NDIMS][VECLEN], double res[VECLEN]) {
    double dotx[VECLEN] ALIGNED;
    dotv(x,x,dotx);
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        res[v] = sqrt(dotx[v]);
    }
}
