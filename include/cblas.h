#ifndef CBLAS_H
#define CBLAS_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


typedef enum CBLAS_LAYOUT {
    CblasRowMajor = 101,
    CblasColMajor = 102
} CBLAS_LAYOUT;

typedef enum CBLAS_TRANSPOSE {
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113,
    CblasConjNoTrans = 114
} CBLAS_TRANSPOSE;

typedef enum CBLAS_UPLO {
    CblasUpper = 121,
    CblasLower = 122
} CBLAS_UPLO;

typedef enum CBLAS_DIAG {
    CblasNonUnit = 131,
    CblasUnit = 132
} CBLAS_DIAG;

typedef enum CBLAS_SIDE {
    CblasLeft = 141,
    CblasRight = 142
} CBLAS_SIDE;


float  cblas_sasum (const int n, const float  *x, const int incx);
double cblas_dasum (const int n, const double *x, const int incx);
float  cblas_scasum(const int n, const void   *x, const int incx);
double cblas_dzasum(const int n, const void   *x, const int incx);

void cblas_saxpy(const int n, const float   alpha, const float  *x, const int incx, float  *y, const int incy);
void cblas_daxpy(const int n, const double  alpha, const double *x, const int incx, double *y, const int incy);
void cblas_caxpy(const int n, const void   *alpha, const void   *x, const int incx, void   *y, const int incy);
void cblas_zaxpy(const int n, const void   *alpha, const void   *x, const int incx, void   *y, const int incy);

void cblas_scopy(const int n, const float  *x, const int incx, float  *y, const int incy);
void cblas_dcopy(const int n, const double *x, const int incx, double *y, const int incy);
void cblas_ccopy(const int n, const void   *x, const int incx, void   *y, const int incy);
void cblas_zcopy(const int n, const void   *x, const int incx, void   *y, const int incy);

float  cblas_sdot(const int n, const float  *x, const int incx, const float  *y, const int incy);
double cblas_ddot(const int n, const double *x, const int incx, const double *y, const int incy);

float  _Complex cblas_cdotc(const int n, const void *x, const int incx, const void *y, const int incy);
double _Complex cblas_zdotc(const int n, const void *x, const int incx, const void *y, const int incy);

void cblas_cdotc_sub(const int n, const void *x, const int incx, const void *y, const int incy, void *ret);
void cblas_zdotc_sub(const int n, const void *x, const int incx, const void *y, const int incy, void *ret);

float  _Complex cblas_cdotu(const int n, const void *x, const int incx, const void *y, const int incy);
double _Complex cblas_zdotu(const int n, const void *x, const int incx, const void *y, const int incy);

void cblas_cdotu_sub(const int n, const void *x, const int incx, const void *y, const int incy, void *ret);
void cblas_zdotu_sub(const int n, const void *x, const int incx, const void *y, const int incy, void *ret);

float  cblas_snrm2 (const int n, const float  *x, const int incx);
double cblas_dnrm2 (const int n, const double *x, const int incx);
float  cblas_scnrm2(const int n, const void   *x, const int incx);
double cblas_dznrm2(const int n, const void  *x, const int incx);

void cblas_srot (const int n, float  *x, const int incx, float  *y, const int incy, const float  c, const float  s);
void cblas_drot (const int n, double *x, const int incx, double *y, const int incy, const double c, const double s);
void cblas_csrot(const int n, void   *x, const int incx, void   *y, const int incy, const float  c, const float  s);
void cblas_zdrot(const int n, void   *x, const int incx, void   *y, const int incy, const double c, const double s);

void cblas_srotg(float  *a, float  *b, float  *c, float  *s);
void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_crotg(void   *a, void   *b, float  *c, void   *s);
void cblas_zrotg(void   *a, void   *b, double *c, void   *s);

void cblas_srotm(const int n, float  *x, const int incx, float  *y, const int incy, const float  *param);
void cblas_drotm(const int n, double *x, const int incx, double *y, const int incy, const double *param);

void cblas_srotmg(float  *d1, float  *d2, float  *x1, const float  y1, float  *param);
void cblas_drotmg(double *d1, double *d2, double *x1, const double y1, double *param);

void cblas_sscal(const int n, const float   alpha, float  *x, const int incx);
void cblas_dscal(const int n, const double  alpha, double *x, const int incx);
void cblas_cscal(const int n, const void   *alpha, void   *x, const int incx);
void cblas_zscal(const int n, const void   *alpha, void   *x, const int incx);
void cblas_csscal(const int n, const float  alpha, void *x, const int incx);
void cblas_zdscal(const int n, const double alpha, void *x, const int incx);

void cblas_sswap(const int n, float  *x, const int incx, float  *y, const int incy);
void cblas_dswap(const int n, double *x, const int incx, double *y, const int incy);
void cblas_cswap(const int n, void   *x, const int incx, void   *y, const int incy);
void cblas_zswap(const int n, void   *x, const int incx, void   *y, const int incy);

unsigned int cblas_isamax(const int n, const float  *x, const int incx);
unsigned int cblas_idamax(const int n, const double *x, const int incx);
unsigned int cblas_icamax(const int n, const void   *x, const int incx);
unsigned int cblas_izamax(const int n, const void   *x, const int incx);

unsigned int cblas_isamin(const int n, const float  *x, const int incx);
unsigned int cblas_idamin(const int n, const double *x, const int incx);
unsigned int cblas_icamin(const int n, const void   *x, const int incx);
unsigned int cblas_izamin(const int n, const void   *x, const int incx);


void cblas_sgemv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const int m, const int n, const float   alpha, const float  *a, const int lda, const float  *x, const int incx, const float   beta, float  *y, const int incy);
void cblas_dgemv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const int m, const int n, const double  alpha, const double *a, const int lda, const double *x, const int incx, const double  beta, double *y, const int incy);
void cblas_cgemv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const int m, const int n, const void   *alpha, const void   *a, const int lda, const void   *x, const int incx, const void   *beta, void   *y, const int incy);
void cblas_zgemv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const int m, const int n, const void   *alpha, const void   *a, const int lda, const void   *x, const int incx, const void   *beta, void   *y, const int incy);

void cblas_sger(const CBLAS_LAYOUT layout, const int m, const int n, const float  alpha, const float  *x, const int incx, const float  *y, const int incy, float  *a, const int lda);
void cblas_dger(const CBLAS_LAYOUT layout, const int m, const int n, const double alpha, const double *x, const int incx, const double *y, const int incy, double *a, const int lda);

void cblas_cgerc(const CBLAS_LAYOUT layout, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *a, const int lda);
void cblas_zgerc(const CBLAS_LAYOUT layout, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *a, const int lda);

void cblas_cgeru(const CBLAS_LAYOUT layout, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *a, const int lda);
void cblas_zgeru(const CBLAS_LAYOUT layout, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *a, const int lda);

void cblas_chemv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *a, const int lda, const void *x, const int incx, const void *beta, void *y, const int incy);
void cblas_zhemv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *a, const int lda, const void *x, const int incx, const void *beta, void *y, const int incy);

void cblas_cher(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const float  alpha, const void *x, const int incx, void *a, const int lda);
void cblas_zher(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const double alpha, const void *x, const int incx, void *a, const int lda);

void cblas_cher2(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *a, const int lda);
void cblas_zher2(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *a, const int lda);

void cblas_ssymv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *a, const int lda, const float  *x, const int incx, const float  beta, float  *y, const int incy);
void cblas_dsymv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const double alpha, const double *a, const int lda, const double *x, const int incx, const double beta, double *y, const int incy);

void cblas_ssyr(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *x, const int incx, float  *a, const int lda);
void cblas_dsyr(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const double alpha, const double *x, const int incx, double *a, const int lda);

void cblas_ssyr2(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *x, const int incx, const float  *y, const int incy, float  *a, const int lda);
void cblas_dsyr2(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const int n, const double alpha, const double *x, const int incx, const double *y, const int incy, double *a, const int lda);

void cblas_strmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const float  *a, const int lda, float  *x, const int incx);
void cblas_dtrmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const double *a, const int lda, double *x, const int incx);
void cblas_ctrmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const void   *a, const int lda, void   *x, const int incx);
void cblas_ztrmv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const void   *a, const int lda, void   *x, const int incx);

void cblas_strsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const float  *a, const int lda, float  *x, const int incx);
void cblas_dtrsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const double *a, const int lda, double *x, const int incx);
void cblas_ctrsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const void   *a, const int lda, void   *x, const int incx);
void cblas_ztrsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int n, const void   *a, const int lda, void   *x, const int incx);


void cblas_sgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k, const float   alpha, const float  *a, const int lda, const float  *b, const int ldb, const float   beta, float  *c, const int ldc);
void cblas_dgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k, const double  alpha, const double *a, const int lda, const double *b, const int ldb, const double  beta, double *c, const int ldc);
void cblas_cgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k, const void   *alpha, const void   *a, const int lda, const void   *b, const int ldb, const void   *beta, void   *c, const int ldc);
void cblas_zgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k, const void   *alpha, const void   *a, const int lda, const void   *b, const int ldb, const void   *beta, void   *c, const int ldc);

void cblas_chemm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void *alpha, const void *a, const int lda, const void *b, const int ldb, const void *beta, void *c, const int ldc);
void cblas_zhemm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void *alpha, const void *a, const int lda, const void *b, const int ldb, const void *beta, void *c, const int ldc);

void cblas_cherk(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const float  alpha, const void *a, const int lda, const float  beta, void *c, const int ldc);
void cblas_zherk(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const double alpha, const void *a, const int lda, const double beta, void *c, const int ldc);

void cblas_cher2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void *alpha, const void *a, const int lda, const void *b, const int ldb, const float  beta, void *c, const int ldc);
void cblas_zher2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void *alpha, const void *a, const int lda, const void *b, const int ldb, const double beta, void *c, const int ldc);

void cblas_ssymm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const float   alpha, const float  *a, const int lda, const float  *b, const int ldb, const float   beta, float   *c, const int ldc);
void cblas_dsymm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const double  alpha, const double *a, const int lda, const double *b, const int ldb, const double  beta, double  *c, const int ldc);
void cblas_csymm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void   *alpha, const void   *a, const int lda, const void   *b, const int ldb, const void   *beta, void    *c, const int ldc);
void cblas_zsymm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void   *alpha, const void   *a, const int lda, const void   *b, const int ldb, const void   *beta, void    *c, const int ldc);

void cblas_ssyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const float   alpha, const float  *a, const int lda, const float   beta, float  *c, const int ldc);
void cblas_dsyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const double  alpha, const double *a, const int lda, const double  beta, double *c, const int ldc);
void cblas_csyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void   *alpha, const void   *a, const int lda, const void   *beta, void   *c, const int ldc);
void cblas_zsyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void   *alpha, const void   *a, const int lda, const void   *beta, void   *c, const int ldc);

void cblas_ssyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const float   alpha, const float  *a, const int lda, const float  *b, const int ldb, const float   beta, float  *c, const int ldc);
void cblas_dsyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const double  alpha, const double *a, const int lda, const double *b, const int ldb, const double  beta, double *c, const int ldc);
void cblas_csyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void   *alpha, const void   *a, const int lda, const void   *b, const int ldb, const void   *beta, void   *c, const int ldc);
void cblas_zsyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void   *alpha, const void   *a, const int lda, const void   *b, const int ldb, const void   *beta, void   *c, const int ldc);

void cblas_strmm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const float   alpha, const float  *a, const int lda, float  *b, const int ldb);
void cblas_dtrmm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const double  alpha, const double *a, const int lda, double *b, const int ldb);
void cblas_ctrmm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const void   *alpha, const void   *a, const int lda, void   *b, const int ldb);
void cblas_ztrmm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const void   *alpha, const void   *a, const int lda, void   *b, const int ldb);

void cblas_strsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const float   alpha, const float  *a, const int lda, float  *b, const int ldb);
void cblas_dtrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const double  alpha, const double *a, const int lda, double *b, const int ldb);
void cblas_ctrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const void   *alpha, const void   *a, const int lda, void   *b, const int ldb);
void cblas_ztrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const int m, const int n, const void   *alpha, const void   *a, const int lda, void   *b, const int ldb);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* CBLAS_H */
