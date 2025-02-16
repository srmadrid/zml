#ifndef CBLAS_H
#define CBLAS_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


typedef enum CBLAS_ORDER {
    CblasRowMajor = 101,
    CblasColMajor = 102
} CBLAS_ORDER;

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

typedef CBLAS_ORDER CBLAS_LAYOUT;


/*
 * Level 1
 */
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


/*
 * Level 2
 */
void cblas_sgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const int kl, const int ku, const float   alpha, const float  *A, const int lda, const float  *x, const int incx, const float   beta, float  *y, const int incy);
void cblas_dgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const int kl, const int ku, const double  alpha, const double *A, const int lda, const double *x, const int incx, const double  beta, double *y, const int incy);
void cblas_cgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const int kl, const int ku, const void   *alpha, const void   *A, const int lda, const void   *x, const int incx, const void   *beta, void   *y, const int incy);
void cblas_zgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const int kl, const int ku, const void   *alpha, const void   *A, const int lda, const void   *x, const int incx, const void   *beta, void   *y, const int incy);

void cblas_sgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const float   alpha, const float  *A, const int lda, const float  *x, const int incx, const float   beta, float  *y, const int incy);
void cblas_dgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const double  alpha, const double *A, const int lda, const double *x, const int incx, const double  beta, double *y, const int incy);
void cblas_cgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const void   *alpha, const void   *A, const int lda, const void   *x, const int incx, const void   *beta, void   *y, const int incy);
void cblas_zgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const int m, const int n, const void   *alpha, const void   *A, const int lda, const void   *x, const int incx, const void   *beta, void   *y, const int incy);

void cblas_sger(const CBLAS_ORDER order, const int m, const int n, const float  alpha, const float  *x, const int incx, const float  *y, const int incy, float  *A, const int lda);
void cblas_dger(const CBLAS_ORDER order, const int m, const int n, const double alpha, const double *x, const int incx, const double *y, const int incy, double *A, const int lda);

void cblas_cgerc(const CBLAS_ORDER order, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *A, const int lda);
void cblas_zgerc(const CBLAS_ORDER order, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *A, const int lda);

void cblas_cgeru(const CBLAS_ORDER order, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *A, const int lda);
void cblas_zgeru(const CBLAS_ORDER order, const int m, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *A, const int lda);

void cblas_chbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const int k, const void *alpha, const void *A, const int lda, const void *x, const int incx, const void *beta, void *y, const int incy);
void cblas_zhbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const int k, const void *alpha, const void *A, const int lda, const void *x, const int incx, const void *beta, void *y, const int incy);

void cblas_chemv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *A, const int lda, const void *x, const int incx, const void *beta, void *y, const int incy);
void cblas_zhemv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *A, const int lda, const void *x, const int incx, const void *beta, void *y, const int incy);

void cblas_cher(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const void *x, const int incx, void *A, const int lda);
void cblas_zher(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const void *x, const int incx, void *A, const int lda);

void cblas_cher2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *A, const int lda);
void cblas_zher2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *A, const int lda);

void cblas_chpmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *Ap, const void *x, const int incx, const void *beta, void *y, const int incy);
void cblas_zhpmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *Ap, const void *x, const int incx, const void *beta, void *y, const int incy);

void cblas_chpr(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const void *x, const int incx, void *Ap);
void cblas_zhpr(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const void *x, const int incx, void *Ap);

void cblas_chpr2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *Ap);
void cblas_zhpr2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const void *alpha, const void *x, const int incx, const void *y, const int incy, void *Ap);

void cblas_ssbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const int k, const float  alpha, const float  *A, const int lda, const float  *x, const int incx, const float  beta, float  *y, const int incy);
void cblas_dsbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const int k, const double alpha, const double *A, const int lda, const double *x, const int incx, const double beta, double *y, const int incy);

void cblas_sspmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *Ap, const float  *x, const int incx, const float  beta, float  *y, const int incy);
void cblas_dspmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const double *Ap, const double *x, const int incx, const double beta, double *y, const int incy);

void cblas_sspr(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *x, const int incx, float  *Ap);
void cblas_dspr(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const double *x, const int incx, double *Ap);

void cblas_sspr2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *x, const int incx, const float  *y, const int incy, float  *Ap);
void cblas_dspr2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const double *x, const int incx, const double *y, const int incy, double *Ap);

void cblas_ssymv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *A, const int lda, const float  *x, const int incx, const float  beta, float  *y, const int incy);
void cblas_dsymv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const double *A, const int lda, const double *x, const int incx, const double beta, double *y, const int incy);

void cblas_ssyr(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *x, const int incx, float  *A, const int lda);
void cblas_dsyr(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const double *x, const int incx, double *A, const int lda);

void cblas_ssyr2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const float  alpha, const float  *x, const int incx, const float  *y, const int incy, float  *A, const int lda);
void cblas_dsyr2(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const int n, const double alpha, const double *x, const int incx, const double *y, const int incy, double *A, const int lda);

void cblas_stbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const float  *A, const int lda, float  *x, const int incx);
void cblas_dtbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const double *A, const int lda, double *x, const int incx);
void cblas_ctbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const void   *A, const int lda, void   *x, const int incx);
void cblas_ztbmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const void   *A, const int lda, void   *x, const int incx);

void cblas_stbsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const float  *A, const int lda, float  *x, const int incx);
void cblas_dtbsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const double *A, const int lda, double *x, const int incx);
void cblas_ctbsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const void   *A, const int lda, void   *x, const int incx);
void cblas_ztbsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const int k, const void   *A, const int lda, void   *x, const int incx);

void cblas_stpmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const float  *Ap, float  *x, const int incx);
void cblas_dtpmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const double *Ap, double *x, const int incx);
void cblas_ctpmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *Ap, void   *x, const int incx);
void cblas_ztpmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *Ap, void   *x, const int incx);

void cblas_stpsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const float  *Ap, float  *x, const int incx);
void cblas_dtpsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const double *Ap, double *x, const int incx);
void cblas_ctpsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *Ap, void   *x, const int incx);
void cblas_ztpsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *Ap, void   *x, const int incx);

void cblas_strmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const float  *A, const int lda, float  *x, const int incx);
void cblas_dtrmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const double *A, const int lda, double *x, const int incx);
void cblas_ctrmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *A, const int lda, void   *x, const int incx);
void cblas_ztrmv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *A, const int lda, void   *x, const int incx);

void cblas_strsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const float  *A, const int lda, float  *x, const int incx);
void cblas_dtrsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const double *A, const int lda, double *x, const int incx);
void cblas_ctrsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *A, const int lda, void   *x, const int incx);
void cblas_ztrsv(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int n, const void   *A, const int lda, void   *x, const int incx);


/*
 * Level 3
 */
void cblas_sgemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, const int m, const int n, const int k, const float   alpha, const float  *A, const int lda, const float  *B, const int ldb, const float   beta, float  *C, const int ldc);
void cblas_dgemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, const int m, const int n, const int k, const double  alpha, const double *A, const int lda, const double *B, const int ldb, const double  beta, double *C, const int ldc);
void cblas_cgemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, const int m, const int n, const int k, const void   *alpha, const void   *A, const int lda, const void   *B, const int ldb, const void   *beta, void   *C, const int ldc);
void cblas_zgemm(const CBLAS_ORDER order, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, const int m, const int n, const int k, const void   *alpha, const void   *A, const int lda, const void   *B, const int ldb, const void   *beta, void   *C, const int ldc);

void cblas_chemm(const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc);
void cblas_zhemm(const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void *alpha, const void *A, const int lda, const void *B, const int ldb, const void *beta, void *C, const int ldc);

void cblas_cherk(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const float  alpha, const void *A, const int lda, const float  beta, void *C, const int ldc);
void cblas_zherk(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const double alpha, const void *A, const int lda, const double beta, void *C, const int ldc);

void cblas_cher2k(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void *alpha, const void *A, const int lda, const void *B, const int ldb, const float  beta, void *C, const int ldc);
void cblas_zher2k(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void *alpha, const void *A, const int lda, const void *B, const int ldb, const double beta, void *C, const int ldc);

void cblas_ssymm(const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const float   alpha, const float  *A, const int lda, const float  *B, const int ldb, const float   beta, float   *C, const int ldc);
void cblas_dsymm(const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const double  alpha, const double *A, const int lda, const double *B, const int ldb, const double  beta, double  *C, const int ldc);
void cblas_csymm(const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void   *alpha, const void   *A, const int lda, const void   *B, const int ldb, const void   *beta, void    *C, const int ldc);
void cblas_zsymm(const CBLAS_ORDER order, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const int m, const int n, const void   *alpha, const void   *A, const int lda, const void   *B, const int ldb, const void   *beta, void    *C, const int ldc);

void cblas_ssyrk(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const float   alpha, const float  *A, const int lda, const float   beta, float  *C, const int ldc);
void cblas_dsyrk(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const double  alpha, const double *A, const int lda, const double  beta, double *C, const int ldc);
void cblas_csyrk(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void   *alpha, const void   *A, const int lda, const void   *beta, void   *C, const int ldc);
void cblas_zsyrk(const CBLAS_ORDER order, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const int n, const int k, const void   *alpha, const void   *A, const int lda, const void   *beta, void   *C, const int ldc);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* CBLAS_H */
