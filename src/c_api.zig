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
    return zml.blas.asum(f32, n, x, incx);
}
export fn cblas_dasum(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.blas.asum(f64, n, x, incx);
}
export fn cblas_casum(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.blas.asum(Complex(f32), n, x, incx);
}
export fn cblas_zasum(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.blas.asum(Complex(f64), n, x, incx);
}

export fn cblas_saxpy(n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.blas.axpy(f32, n, alpha, x, incx, y, incy);
}
export fn cblas_daxpy(n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.blas.axpy(f64, n, alpha, x, incx, y, incy);
}
export fn cblas_caxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.axpy(Complex(f32), n, alpha.*, x, incx, y, incy);
}
export fn cblas_zaxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.axpy(Complex(f64), n, alpha.*, x, incx, y, incy);
}

export fn cblas_scopy(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.blas.copy(f32, n, x, incx, y, incy);
}
export fn cblas_dcopy(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.blas.copy(f64, n, x, incx, y, incy);
}
export fn cblas_ccopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.copy(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zcopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.copy(Complex(f64), n, x, incx, y, incy);
}

export fn cblas_sdot(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int) f32 {
    return zml.blas.dot(f32, n, x, incx, y, incy);
}
export fn cblas_ddot(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int) f64 {
    return zml.blas.dot(f64, n, x, incx, y, incy);
}

export fn cblas_cdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f32) {
    return zml.blas.dotc(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f64) {
    return zml.blas.dotc(Complex(f64), n, x, incx, y, incy);
}
export fn cblas_cdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.blas.dotc_sub(Complex(f32), n, x, incx, y, incy, ret);
}
export fn cblas_zdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.blas.dotc_sub(Complex(f64), n, x, incx, y, incy, ret);
}

export fn cblas_cdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f32) {
    return zml.blas.dotu(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) Complex(f64) {
    return zml.blas.dotu(Complex(f64), n, x, incx, y, incy);
}
export fn cblas_cdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.blas.dotu_sub(Complex(f32), n, x, incx, y, incy, ret);
}
export fn cblas_zdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.blas.dotu_sub(Complex(f64), n, x, incx, y, incy, ret);
}

export fn cblas_snrm2(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.blas.nrm2(f32, n, x, incx);
}
export fn cblas_dnrm2(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.blas.nrm2(f64, n, x, incx);
}
export fn cblas_scnrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.blas.nrm2(Complex(f32), n, x, incx);
}
export fn cblas_dznrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.blas.nrm2(Complex(f64), n, x, incx);
}

export fn cblas_srot(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, c: f32, s: f32) void {
    return zml.blas.rot(f32, n, x, incx, y, incy, c, s);
}
export fn cblas_drot(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, c: f64, s: f64) void {
    return zml.blas.rot(f64, n, x, incx, y, incy, c, s);
}
export fn cblas_csrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f32, s: f32) void {
    return zml.blas.rot(Complex(f32), n, x, incx, y, incy, c, s);
}
export fn cblas_zdrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f64, s: f64) void {
    return zml.blas.rot(Complex(f64), n, x, incx, y, incy, c, s);
}

export fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return zml.blas.rotg(f32, a, b, c, s);
}
export fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return zml.blas.rotg(f64, a, b, c, s);
}
export fn cblas_crotg(a: *anyopaque, b: *anyopaque, c: *f32, s: *anyopaque) void {
    return zml.blas.rotg(Complex(f32), a, b, c, s);
}
export fn cblas_zrotg(a: *anyopaque, b: *anyopaque, c: *f64, s: *anyopaque) void {
    return zml.blas.rotg(Complex(f64), a, b, c, s);
}

export fn cblas_srotm(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, param: [*c]f32) void {
    return zml.blas.rotm(f32, n, x, incx, y, incy, param);
}
export fn cblas_drotm(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, param: [*c]f64) void {
    return zml.blas.rotm(f64, n, x, incx, y, incy, param);
}

export fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void {
    return zml.blas.rotmg(f32, d1, d2, x1, y1, param);
}
export fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void {
    return zml.blas.rotmg(f64, d1, d2, x1, y1, param);
}

export fn cblas_sscal(n: c_int, alpha: f32, x: [*c]f32, incx: c_int) void {
    return zml.blas.scal(f32, n, alpha, x, incx);
}
export fn cblas_dscal(n: c_int, alpha: f64, x: [*c]f64, incx: c_int) void {
    return zml.blas.scal(f64, n, alpha, x, incx);
}
export fn cblas_cscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.scal(Complex(f32), n, alpha.*, x, incx);
}
export fn cblas_zscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.scal(Complex(f64), n, alpha.*, x, incx);
}
export fn cblas_csscal(n: c_int, alpha: f32, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.csscal(n, alpha, x, incx);
}
export fn cblas_zdscal(n: c_int, alpha: f64, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.zdscal(n, alpha, x, incx);
}

export fn cblas_sswap(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.blas.swap(f32, n, x, incx, y, incy);
}
export fn cblas_dswap(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.blas.swap(f64, n, x, incx, y, incy);
}
export fn cblas_cswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.swap(Complex(f32), n, x, incx, y, incy);
}
export fn cblas_zswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.swap(Complex(f64), n, x, incx, y, incy);
}

export fn cblas_isamax(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.blas.iamax(f32, n, x, incx);
}
export fn cblas_idamax(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.blas.iamax(f64, n, x, incx);
}
export fn cblas_icamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.blas.iamax(Complex(f32), n, x, incx);
}
export fn cblas_izamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.blas.iamax(Complex(f64), n, x, incx);
}

export fn cblas_isamin(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.blas.iamin(f32, n, x, incx);
}
export fn cblas_idamin(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.blas.iamin(f64, n, x, incx);
}
export fn cblas_icamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.blas.iamin(Complex(f32), n, x, incx);
}
export fn cblas_izamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.blas.iamin(Complex(f64), n, x, incx);
}

export fn cblas_sgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.blas.gbmv(f32, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.blas.gbmv(f64, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_cgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.gbmv(Complex(f32), order, transA, m, n, kl, ku, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zgbmv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.gbmv(Complex(f64), order, transA, m, n, kl, ku, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_sgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.blas.gemv(f32, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.blas.gemv(f64, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_cgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.gemv(Complex(f32), order, transA, m, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zgemv(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.gemv(Complex(f64), order, transA, m, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_sger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, A: [*c]f32, lda: c_int) void {
    return zml.blas.ger(f32, order, m, n, alpha, x, incx, y, incy, A, lda);
}
export fn cblas_dger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, A: [*c]f64, lda: c_int) void {
    return zml.blas.ger(f64, order, m, n, alpha, x, incx, y, incy, A, lda);
}

export fn cblas_cgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.geru(Complex(f32), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}
export fn cblas_zgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.geru(Complex(f64), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}

export fn cblas_cgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.gerc(Complex(f32), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}
export fn cblas_zgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.gerc(Complex(f64), order, m, n, alpha.*, x, incx, y, incy, A, lda);
}

export fn cblas_chbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.hbmv(Complex(f32), order, uplo, n, k, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zhbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.hbmv(Complex(f64), order, uplo, n, k, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_chemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.hemv(Complex(f32), order, uplo, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}
export fn cblas_zhemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.hemv(Complex(f64), order, uplo, n, alpha.*, A, lda, x, incx, beta.*, y, incy);
}

export fn cblas_cher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.her(Complex(f32), order, uplo, n, alpha, x, incx, A, lda);
}
export fn cblas_zher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.her(Complex(f64), order, uplo, n, alpha, x, incx, A, lda);
}

export fn cblas_cher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.her2(Complex(f32), order, uplo, n, alpha.*, x, incx, y, incy, A, lda);
}
export fn cblas_zher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, A: [*c]anyopaque, lda: c_int) void {
    return zml.blas.her2(Complex(f64), order, uplo, n, alpha.*, x, incx, y, incy, A, lda);
}

export fn cblas_chpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.hpmv(Complex(f32), order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy);
}
export fn cblas_zhpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.blas.hpmv(Complex(f64), order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy);
}

export fn cblas_chpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.blas.hpr(Complex(f32), order, uplo, n, alpha, x, incx, Ap);
}
export fn cblas_zhpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.blas.hpr(Complex(f64), order, uplo, n, alpha, x, incx, Ap);
}

export fn cblas_chpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.blas.hpr2(Complex(f32), order, uplo, n, alpha.*, x, incx, y, incy, Ap);
}
export fn cblas_zhpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.blas.hpr2(Complex(f64), order, uplo, n, alpha.*, x, incx, y, incy, Ap);
}

export fn cblas_ssbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.blas.sbmv(f32, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dsbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.blas.sbmv(f64, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

export fn cblas_sspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, Ap: [*c]const f32, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.blas.spmv(f32, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
export fn cblas_dspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, Ap: [*c]const f64, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.blas.spmv(f64, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}

export fn cblas_sspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, Ap: [*c]f32) void {
    return zml.blas.spr(f32, order, uplo, n, alpha, x, incx, Ap);
}
export fn cblas_dspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, Ap: [*c]f64) void {
    return zml.blas.spr(f64, order, uplo, n, alpha, x, incx, Ap);
}

export fn cblas_sspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, Ap: [*c]f32) void {
    return zml.blas.spr2(f32, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
export fn cblas_dspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, Ap: [*c]f64) void {
    return zml.blas.spr2(f64, order, uplo, n, alpha, x, incx, y, incy, Ap);
}

export fn cblas_ssymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, A: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.blas.symv(f32, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
export fn cblas_dsymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, A: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.blas.symv(f64, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}

export fn cblas_ssyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, A: [*c]f32, lda: c_int) void {
    return zml.blas.syr(f32, order, uplo, n, alpha, x, incx, A, lda);
}
export fn cblas_dsyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, A: [*c]f64, lda: c_int) void {
    return zml.blas.syr(f64, order, uplo, n, alpha, x, incx, A, lda);
}

export fn cblas_ssyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, A: [*c]f32, lda: c_int) void {
    return zml.blas.syr2(f32, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
export fn cblas_dsyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, A: [*c]f64, lda: c_int) void {
    return zml.blas.syr2(f64, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}

export fn cblas_stbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.blas.tbmv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_dtbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.blas.tbmv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ctbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tbmv(Complex(f32), order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ztbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tbmv(Complex(f64), order, uplo, transA, diag, n, k, A, lda, x, incx);
}

export fn cblas_stbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.blas.tbsv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_dtbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.blas.tbsv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ctbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tbsv(Complex(f32), order, uplo, transA, diag, n, k, A, lda, x, incx);
}
export fn cblas_ztbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tbsv(Complex(f64), order, uplo, transA, diag, n, k, A, lda, x, incx);
}

export fn cblas_stpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.blas.tpmv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_dtpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.blas.tpmv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_ctpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tpmv(Complex(f32), order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_ztpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tpmv(Complex(f64), order, uplo, transA, diag, n, Ap, x, incx);
}

export fn cblas_stpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.blas.tpsv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_dtpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.blas.tpsv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_ctpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tpsv(Complex(f32), order, uplo, transA, diag, n, Ap, x, incx);
}
export fn cblas_ztpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.tpsv(Complex(f64), order, uplo, transA, diag, n, Ap, x, incx);
}

export fn cblas_strmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.blas.trmv(f32, order, uplo, transA, diag, n, A, lda, x, incx);
}
export fn cblas_dtrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.blas.trmv(f64, order, uplo, transA, diag, n, A, lda, x, incx);
}
export fn cblas_ctrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.trmv(Complex(f32), order, uplo, transA, diag, n, A, lda, x, incx);
}
export fn cblas_ztrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.trmv(Complex(f64), order, uplo, transA, diag, n, A, lda, x, incx);
}

export fn cblas_strsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.blas.trsv(f32, order, uplo, transA, diag, n, A, lda, x, incx);
}
export fn cblas_dtrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.blas.trsv(f64, order, uplo, transA, diag, n, A, lda, x, incx);
}
export fn cblas_ctrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.trsv(Complex(f32), order, uplo, transA, diag, n, A, lda, x, incx);
}
export fn cblas_ztrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transA: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, A: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.blas.trsv(Complex(f64), order, uplo, transA, diag, n, A, lda, x, incx);
}

export fn cblas_sgemm(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, transB: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f32, A: [*c]const f32, lda: c_int, B: [*c]const f32, ldb: c_int, beta: f32, C: [*c]f32, ldc: c_int) void {
    return zml.blas.gemm(f32, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
export fn cblas_dgemm(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, transB: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f64, A: [*c]const f64, lda: c_int, B: [*c]const f64, ldb: c_int, beta: f64, C: [*c]f64, ldc: c_int) void {
    return zml.blas.gemm(f64, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
export fn cblas_cgemm(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, transB: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, B: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.gemm(Complex(f32), order, transA, transB, m, n, k, alpha.*, A, lda, B, ldb, beta.*, C, ldc);
}
export fn cblas_zgemm(order: CBLAS_ORDER, transA: CBLAS_TRANSPOSE, transB: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, B: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.gemm(Complex(f64), order, transA, transB, m, n, k, alpha.*, A, lda, B, ldb, beta.*, C, ldc);
}

export fn cblas_chemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, B: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.hemm(Complex(f32), order, side, uplo, m, n, alpha.*, A, lda, B, ldb, beta.*, C, ldc);
}
export fn cblas_zhemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, B: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.hemm(Complex(f64), order, side, uplo, m, n, alpha.*, A, lda, B, ldb, beta.*, C, ldc);
}

export fn cblas_cherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, A: [*c]const anyopaque, lda: c_int, beta: f32, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.herk(Complex(f32), order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
export fn cblas_zherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, A: [*c]const anyopaque, lda: c_int, beta: f64, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.herk(Complex(f64), order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

export fn cblas_cher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, B: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.her2k(Complex(f32), order, uplo, trans, n, k, alpha.*, A, lda, B, ldb, beta.*, C, ldc);
}
export fn cblas_zher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, A: [*c]const anyopaque, lda: c_int, B: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, C: [*c]anyopaque, ldc: c_int) void {
    return zml.blas.her2k(Complex(f64), order, uplo, trans, n, k, alpha.*, A, lda, B, ldb, beta.*, C, ldc);
}
