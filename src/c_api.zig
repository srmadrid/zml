const std = @import("std");
const zml = @import("zml");
const Complex = std.math.Complex;

// CBLAS
export const CBLAS_ORDER = enum(c_int) {
    CblasRowMajor = 101,
    CblasColMajor = 102,
};

export const CBLAS_TRANSPOSE = enum(c_int) {
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113,
    CblasConjNoTrans = 114,
};

export const CBLAS_UPLO = enum(c_int) {
    CblasUpper = 121,
    CblasLower = 122,
};

export const CBLAS_DIAG = enum(c_int) {
    CblasNonUnit = 131,
    CblasUnit = 132,
};

export const CBLAS_SIDE = enum(c_int) {
    CblasLeft = 141,
    CblasRight = 142,
};

export fn cblas_sasum(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.BLAS.asum(f32, n, x, incx);
}
export fn cblas_dasum(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.BLAS.asum(f64, n, x, incx);
}
export fn cblas_casum(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.BLAS.asum(Complex(f32), n, x, incx);
}
export fn cblas_zasum(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.BLAS.asum(Complex(f64), n, x, incx);
}

export fn cblas_saxpy(n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.axpy(f32, n, alpha, x, incx, y, incy);
}
export fn cblas_daxpy(n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.axpy(f64, n, alpha, x, incx, y, incy);
}
export fn cblas_caxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.axpy(Complex(f32), n, alpha.*, x, incx, y, incy);
}
export fn cblas_zaxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.axpy(Complex(f64), n, alpha.*, x, incx, y, incy);
}

export fn cblas_scopy(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.copy(f32, n, x, incx, y, incy);
}
export fn cblas_dcopy(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.copy(f64, n, x, incx, y, incy);
}
export fn cblas_ccopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.copy(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zcopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.copy(Complex(f64), n, x, incx, y, incy);
}

export fn cblas_sdot(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int) f32 {
    return zml.BLAS.dot(f32, n, x, incx, y, incy);
}
export fn cblas_ddot(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int) f64 {
    return zml.BLAS.dot(f64, n, x, incx, y, incy);
}

export fn cblas_cdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f32) {
    return zml.BLAS.dotc(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f64) {
    return zml.BLAS.dotc(Complex(f64), n, x, incx, y, incy);
}
export fn cblas_cdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.BLAS.dotc_sub(Complex(f32), n, x, incx, y, incy, ret);
}
export fn cblas_zdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.BLAS.dotc_sub(Complex(f64), n, x, incx, y, incy, ret);
}

export fn cblas_cdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f32) {
    return zml.BLAS.dotu(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f64) {
    return zml.BLAS.dotu(Complex(f64), n, x, incx, y, incy);
}
export fn cblas_cdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.BLAS.dotu_sub(Complex(f32), n, x, incx, y, incy, ret);
}
export fn cblas_zdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.BLAS.dotu_sub(Complex(f64), n, x, incx, y, incy, ret);
}

export fn cblas_snrm2(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.BLAS.nrm2(f32, n, x, incx);
}
export fn cblas_dnrm2(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.BLAS.nrm2(f64, n, x, incx);
}
export fn cblas_scnrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.BLAS.nrm2(Complex(f32), n, x, incx);
}
export fn cblas_dznrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.BLAS.nrm2(Complex(f64), n, x, incx);
}

export fn cblas_srot(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, c: f32, s: f32) void {
    return zml.BLAS.rot(f32, n, x, incx, y, incy, c, s);
}
export fn cblas_drot(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, c: f64, s: f64) void {
    return zml.BLAS.rot(f64, n, x, incx, y, incy, c, s);
}
export fn cblas_csrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f32, s: f32) void {
    return zml.BLAS.rot(Complex(f32), n, x, incx, y, incy, c, s);
}
export fn cblas_zdrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f64, s: f64) void {
    return zml.BLAS.rot(Complex(f64), n, x, incx, y, incy, c, s);
}

export fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return zml.BLAS.rotg(f32, a, b, c, s);
}
export fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return zml.BLAS.rotg(f64, a, b, c, s);
}
export fn cblas_crotg(a: *anyopaque, b: *anyopaque, c: *f32, s: *anyopaque) void {
    return zml.BLAS.rotg(Complex(f32), a, b, c, s);
}
export fn cblas_zrotg(a: *anyopaque, b: *anyopaque, c: *f64, s: *anyopaque) void {
    return zml.BLAS.rotg(Complex(f64), a, b, c, s);
}

export fn cblas_srotm(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, param: [*c]f32) void {
    return zml.BLAS.rotm(f32, n, x, incx, y, incy, param);
}
export fn cblas_drotm(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, param: [*c]f64) void {
    return zml.BLAS.rotm(f64, n, x, incx, y, incy, param);
}

export fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void {
    return zml.BLAS.rotmg(f32, d1, d2, x1, y1, param);
}
export fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void {
    return zml.BLAS.rotmg(f64, d1, d2, x1, y1, param);
}

export fn cblas_sscal(n: c_int, alpha: f32, x: [*c]f32, incx: c_int) void {
    return zml.BLAS.scal(f32, n, alpha, x, incx);
}
export fn cblas_dscal(n: c_int, alpha: f64, x: [*c]f64, incx: c_int) void {
    return zml.BLAS.scal(f64, n, alpha, x, incx);
}
export fn cblas_cscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.scal(Complex(f32), n, alpha.*, x, incx);
}
export fn cblas_zscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.scal(Complex(f64), n, alpha.*, x, incx);
}

export fn cblas_sswap(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.swap(f32, n, x, incx, y, incy);
}
export fn cblas_dswap(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.swap(f64, n, x, incx, y, incy);
}
export fn cblas_cswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.swap(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.swap(Complex(f64), n, x, incx, y, incy);
}

export fn cblas_isamax(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.BLAS.iamax(f32, n, x, incx);
}
export fn cblas_idamax(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.BLAS.iamax(f64, n, x, incx);
}
export fn cblas_icamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.BLAS.iamax(Complex(f32), n, x, incx);
}
export fn cblas_izamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.BLAS.iamax(Complex(f64), n, x, incx);
}

export fn cblas_isamin(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.BLAS.iamin(f32, n, x, incx);
}
export fn cblas_idamin(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.BLAS.iamin(f64, n, x, incx);
}
export fn cblas_icamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.BLAS.iamin(Complex(f32), n, x, incx);
}
export fn cblas_izamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.BLAS.iamin(Complex(f64), n, x, incx);
}

export fn cblas_sgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.gbmv(f32, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.gbmv(f64, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_cgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.gbmv(Complex(f32), order, transA, m, n, kl, ku, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.gbmv(Complex(f64), order, transA, m, n, kl, ku, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_sgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.gemv(f32, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.gemv(f64, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_cgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.gemv(Complex(f32), order, transA, m, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.gemv(Complex(f64), order, transA, m, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_sger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, A: [*c]f32, lda: c_int) void {
    return zml.BLAS.ger(f32, order, m, n, alpha, x, incx, y, incy, A, lda);
}
export fn cblas_dger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, A: [*c]f64, lda: c_int) void {
    return zml.BLAS.ger(f64, order, m, n, alpha, x, incx, y, incy, A, lda);
}

export fn cblas_cgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.geru(Complex(f32), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}
export fn cblas_zgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.geru(Complex(f64), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}

export fn cblas_cgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.gerc(Complex(f32), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}
export fn cblas_zgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.gerc(Complex(f64), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}

export fn cblas_chbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.hbmv(Complex(f32), order, uplo, n, k, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zhbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.hbmv(Complex(f64), order, uplo, n, k, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_chemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.hemv(Complex(f32), order, uplo, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zhemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.hemv(Complex(f64), order, uplo, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_cher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.her(Complex(f32), order, uplo, n, alpha, x, incx, A, lda);
}
export fn cblas_zher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.her(Complex(f64), order, uplo, n, alpha, x, incx, A, lda);
}

export fn cblas_cher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.her2(Complex(f32), order, uplo, n, alpha.*, x, incx, y, incy, A, lda);
}
export fn cblas_zher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.BLAS.her2(Complex(f64), order, uplo, n, alpha.*, x, incx, y, incy, A, lda);
}

export fn cblas_chpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.hpmv(Complex(f32), order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy);
}
export fn cblas_zhpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.BLAS.hpmv(Complex(f64), order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy);
}

export fn cblas_chpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.BLAS.hpr(Complex(f32), order, uplo, n, alpha, x, incx, Ap);
}
export fn cblas_zhpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.BLAS.hpr(Complex(f64), order, uplo, n, alpha, x, incx, Ap);
}

export fn cblas_chpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.BLAS.hpr2(Complex(f32), order, uplo, n, alpha.*, x, incx, y, incy, Ap);
}
export fn cblas_zhpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.BLAS.hpr2(Complex(f64), order, uplo, n, alpha.*, x, incx, y, incy, Ap);
}

export fn cblas_ssbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.sbmv(f32, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dsbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.sbmv(f64, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

export fn cblas_sspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, Ap: [*c]const f32, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.spmv(f32, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
export fn cblas_dspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, Ap: [*c]const f64, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.spmv(f64, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}

export fn cblas_sspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, Ap: [*c]f32) void {
    return zml.BLAS.spr(f32, order, uplo, n, alpha, x, incx, Ap);
}
export fn cblas_dspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, Ap: [*c]f64) void {
    return zml.BLAS.spr(f64, order, uplo, n, alpha, x, incx, Ap);
}

export fn cblas_sspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, Ap: [*c]f32) void {
    return zml.BLAS.spr2(f32, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
export fn cblas_dspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, Ap: [*c]f64) void {
    return zml.BLAS.spr2(f64, order, uplo, n, alpha, x, incx, y, incy, Ap);
}

export fn cblas_ssymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.BLAS.symv(f32, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dsymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.BLAS.symv(f64, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}

export fn cblas_ssyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, A: [*c]f32, lda: c_int) void {
    return zml.BLAS.syr(f32, order, uplo, n, alpha, x, incx, A, lda);
}
export fn cblas_dsyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, A: [*c]f64, lda: c_int) void {
    return zml.BLAS.syr(f64, order, uplo, n, alpha, x, incx, A, lda);
}

export fn cblas_ssyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, A: [*c]f32, lda: c_int) void {
    return zml.BLAS.syr2(f32, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
export fn cblas_dsyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, A: [*c]f64, lda: c_int) void {
    return zml.BLAS.syr2(f64, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}

export fn cblas_stbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.BLAS.tbmv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_dtbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.BLAS.tbmv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ctbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.tbmv(Complex(f32), order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ztbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.tbmv(Complex(f64), order, uplo, transA, diag, n, k, A, lda, x, incx);
}

export fn cblas_stbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.BLAS.tbsv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_dtbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.BLAS.tbsv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ctbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.tbsv(Complex(f32), order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ztbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.tbsv(Complex(f64), order, uplo, transA, diag, n, k, A, lda, x, incx);
}

export fn cblas_stpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.BLAS.tpmv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_dtpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.BLAS.tpmv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_ctpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.tpmv(Complex(f32), order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_ztpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.BLAS.tpmv(Complex(f64), order, uplo, transA, diag, n, Ap, x, incx);
}
