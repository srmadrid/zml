const std = @import("std");
const zml = @import("zml");

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

// Level 1
export fn cblas_sasum(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.linalg.blas.asum(n, x, incx, .{}) catch unreachable;
}
export fn cblas_dasum(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.linalg.blas.asum(n, x, incx, .{}) catch unreachable;
}
export fn cblas_casum(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.linalg.blas.asum(n, x, incx, .{}) catch unreachable;
}
export fn cblas_zasum(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.linalg.blas.asum(n, x, incx, .{}) catch unreachable;
}

export fn cblas_saxpy(n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.axpy(n, alpha, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_daxpy(n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.axpy(n, alpha, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_caxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.axpy(n, alpha.*, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zaxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.axpy(n, alpha.*, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_scopy(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.copy(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_dcopy(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.copy(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_ccopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.copy(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zcopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.copy(n, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_sdot(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int) f32 {
    return zml.linalg.blas.dot(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_ddot(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int) f64 {
    return zml.linalg.blas.dot(n, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_cdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf32 {
    return zml.linalg.blas.dotc(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf64 {
    return zml.linalg.blas.dotc(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_cdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotc_sub(n, x, incx, y, incy, ret, .{}) catch unreachable;
}
export fn cblas_zdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotc_sub(n, x, incx, y, incy, ret, .{}) catch unreachable;
}

export fn cblas_cdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf32 {
    return zml.linalg.blas.dotu(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf64 {
    return zml.linalg.blas.dotu(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_cdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotu_sub(n, x, incx, y, incy, ret, .{}) catch unreachable;
}
export fn cblas_zdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotu_sub(n, x, incx, y, incy, ret, .{}) catch unreachable;
}

export fn cblas_snrm2(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.linalg.blas.nrm2(n, x, incx, .{}) catch unreachable;
}
export fn cblas_dnrm2(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.linalg.blas.nrm2(n, x, incx, .{}) catch unreachable;
}
export fn cblas_scnrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.linalg.blas.nrm2(n, x, incx, .{}) catch unreachable;
}
export fn cblas_dznrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.linalg.blas.nrm2(n, x, incx, .{}) catch unreachable;
}

export fn cblas_srot(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, c: f32, s: f32) void {
    return zml.linalg.blas.rot(n, x, incx, y, incy, c, s, .{}) catch unreachable;
}
export fn cblas_drot(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, c: f64, s: f64) void {
    return zml.linalg.blas.rot(n, x, incx, y, incy, c, s, .{}) catch unreachable;
}
export fn cblas_csrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f32, s: f32) void {
    return zml.linalg.blas.rot(n, x, incx, y, incy, c, s, .{}) catch unreachable;
}
export fn cblas_zdrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f64, s: f64) void {
    return zml.linalg.blas.rot(n, x, incx, y, incy, c, s, .{}) catch unreachable;
}

export fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return zml.linalg.blas.rotg(a, b, c, s, .{}) catch unreachable;
}
export fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return zml.linalg.blas.rotg(a, b, c, s, .{}) catch unreachable;
}
export fn cblas_crotg(a: *anyopaque, b: *anyopaque, c: *f32, s: *anyopaque) void {
    return zml.linalg.blas.rotg(a, b, c, s, .{}) catch unreachable;
}
export fn cblas_zrotg(a: *anyopaque, b: *anyopaque, c: *f64, s: *anyopaque) void {
    return zml.linalg.blas.rotg(a, b, c, s, .{}) catch unreachable;
}

export fn cblas_srotm(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, param: [*c]f32) void {
    return zml.linalg.blas.rotm(n, x, incx, y, incy, param, .{}) catch unreachable;
}
export fn cblas_drotm(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, param: [*c]f64) void {
    return zml.linalg.blas.rotm(n, x, incx, y, incy, param, .{}) catch unreachable;
}

export fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void {
    return zml.linalg.blas.rotmg(d1, d2, x1, y1, param, .{}) catch unreachable;
}
export fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void {
    return zml.linalg.blas.rotmg(d1, d2, x1, y1, param, .{}) catch unreachable;
}

export fn cblas_sscal(n: c_int, alpha: f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.scal(n, alpha, x, incx, .{}) catch unreachable;
}
export fn cblas_dscal(n: c_int, alpha: f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.scal(n, alpha, x, incx, .{}) catch unreachable;
}
export fn cblas_cscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.scal(n, alpha.*, x, incx, .{}) catch unreachable;
}
export fn cblas_zscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.scal(n, alpha.*, x, incx, .{}) catch unreachable;
}
export fn cblas_csscal(n: c_int, alpha: f32, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.csscal(n, alpha, x, incx, .{}) catch unreachable;
}
export fn cblas_zdscal(n: c_int, alpha: f64, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.zdscal(n, alpha, x, incx, .{}) catch unreachable;
}

export fn cblas_sswap(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.swap(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_dswap(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.swap(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_cswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.swap(n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.swap(n, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_isamax(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(n, x, incx, .{}) catch unreachable;
}
export fn cblas_idamax(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(n, x, incx, .{}) catch unreachable;
}
export fn cblas_icamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(n, x, incx, .{}) catch unreachable;
}
export fn cblas_izamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(n, x, incx, .{}) catch unreachable;
}

export fn cblas_isamin(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(n, x, incx, .{}) catch unreachable;
}
export fn cblas_idamin(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(n, x, incx, .{}) catch unreachable;
}
export fn cblas_icamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(n, x, incx, .{}) catch unreachable;
}
export fn cblas_izamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(n, x, incx, .{}) catch unreachable;
}

// Level 2
export fn cblas_sgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.gbmv(order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_cgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gbmv(order, transa, m, n, kl, ku, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gbmv(order, transa, m, n, kl, ku, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_sgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.gemv(order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_cgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gemv(order, transa, m, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gemv(order, transa, m, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_sger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.ger(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_dger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.ger(order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_cgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.geru(order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_zgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.geru(order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_cgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.gerc(order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_zgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.gerc(order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_chbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hbmv(order, uplo, n, k, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zhbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hbmv(order, uplo, n, k, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_chemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hemv(order, uplo, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zhemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hemv(order, uplo, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_cher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her(order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}
export fn cblas_zher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her(order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}

export fn cblas_cher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her2(order, uplo, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_zher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her2(order, uplo, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_chpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hpmv(order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zhpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hpmv(order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_chpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr(order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}
export fn cblas_zhpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr(order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}

export fn cblas_chpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr2(order, uplo, n, alpha.*, x, incx, y, incy, Ap, .{}) catch unreachable;
}
export fn cblas_zhpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr2(order, uplo, n, alpha.*, x, incx, y, incy, Ap, .{}) catch unreachable;
}

export fn cblas_ssbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.sbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dsbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.sbmv(order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}

export fn cblas_sspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, Ap: [*c]const f32, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.spmv(order, uplo, n, alpha, Ap, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, Ap: [*c]const f64, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.spmv(order, uplo, n, alpha, Ap, x, incx, beta, y, incy, .{}) catch unreachable;
}

export fn cblas_sspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, Ap: [*c]f32) void {
    return zml.linalg.blas.spr(order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}
export fn cblas_dspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, Ap: [*c]f64) void {
    return zml.linalg.blas.spr(order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}

export fn cblas_sspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, Ap: [*c]f32) void {
    return zml.linalg.blas.spr2(order, uplo, n, alpha, x, incx, y, incy, Ap, .{}) catch unreachable;
}
export fn cblas_dspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, Ap: [*c]f64) void {
    return zml.linalg.blas.spr2(order, uplo, n, alpha, x, incx, y, incy, Ap, .{}) catch unreachable;
}

export fn cblas_ssymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.symv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dsymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.symv(order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}

export fn cblas_ssyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.syr(order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}
export fn cblas_dsyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.syr(order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}

export fn cblas_ssyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.syr2(order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_dsyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.syr2(order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_stbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbmv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_stbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbsv(order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_stpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tpmv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_dtpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tpmv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ctpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpmv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ztpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpmv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}

export fn cblas_stpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tpsv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_dtpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tpsv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ctpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpsv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ztpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpsv(order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}

export fn cblas_strmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trmv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_strsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trsv(order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}

// Level 3
export fn cblas_sgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_dgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.gemm(order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_cgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.gemm(order, transa, transb, m, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.gemm(order, transa, transb, m, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_chemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.hemm(order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zhemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.hemm(order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_cherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: [*c]const anyopaque, lda: c_int, beta: f32, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.herk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_zherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: [*c]const anyopaque, lda: c_int, beta: f64, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.herk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}

export fn cblas_cher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.her2k(order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.her2k(order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_ssymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_dsymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.symm(order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_csymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.symm(order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zsymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.symm(order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_ssyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_dsyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.syrk(order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_csyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syrk(order, uplo, trans, n, k, alpha.*, a, lda, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zsyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syrk(order, uplo, trans, n, k, alpha.*, a, lda, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_ssyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_dsyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_csyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zsyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_strmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]f32, ldb: c_int) void {
    return zml.linalg.blas.trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_dtrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]f64, ldb: c_int) void {
    return zml.linalg.blas.trmm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ctrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trmm(order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ztrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trmm(order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}

export fn cblas_strsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]f32, ldb: c_int) void {
    return zml.linalg.blas.trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_dtrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]f64, ldb: c_int) void {
    return zml.linalg.blas.trsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ctrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trsm(order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ztrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trsm(order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}
