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

export fn cblas_sasum(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.linalg.blas.asum(f32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_dasum(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.linalg.blas.asum(f64, n, x, incx, .{}) catch unreachable;
}
export fn cblas_casum(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.linalg.blas.asum(zml.cf32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_zasum(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.linalg.blas.asum(zml.cf64, n, x, incx, .{}) catch unreachable;
}

export fn cblas_saxpy(n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.axpy(f32, n, alpha, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_daxpy(n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.axpy(f64, n, alpha, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_caxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.axpy(zml.cf32, n, alpha.*, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zaxpy(n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.axpy(zml.cf64, n, alpha.*, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_scopy(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.copy(f32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_dcopy(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.copy(f64, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_ccopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.copy(zml.cf32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zcopy(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.copy(zml.cf64, n, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_sdot(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int) f32 {
    return zml.linalg.blas.dot(f32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_ddot(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int) f64 {
    return zml.linalg.blas.dot(f64, n, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_cdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf32 {
    return zml.linalg.blas.dotc(zml.cf32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zdotc(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf64 {
    return zml.linalg.blas.dotc(zml.cf64, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_cdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotc_sub(zml.cf32, n, x, incx, y, incy, ret, .{}) catch unreachable;
}
export fn cblas_zdotc_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotc_sub(zml.cf64, n, x, incx, y, incy, ret, .{}) catch unreachable;
}

export fn cblas_cdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf32 {
    return zml.linalg.blas.dotu(zml.cf32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zdotu(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int) zml.cf64 {
    return zml.linalg.blas.dotu(zml.cf64, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_cdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotu_sub(zml.cf32, n, x, incx, y, incy, ret, .{}) catch unreachable;
}
export fn cblas_zdotu_sub(n: c_int, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.dotu_sub(zml.cf64, n, x, incx, y, incy, ret, .{}) catch unreachable;
}

export fn cblas_snrm2(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.linalg.blas.nrm2(f32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_dnrm2(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.linalg.blas.nrm2(f64, n, x, incx, .{}) catch unreachable;
}
export fn cblas_scnrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f32 {
    return zml.linalg.blas.nrm2(zml.cf32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_dznrm2(n: c_int, x: [*c]const anyopaque, incx: c_int) f64 {
    return zml.linalg.blas.nrm2(zml.cf64, n, x, incx, .{}) catch unreachable;
}

export fn cblas_srot(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, c: f32, s: f32) void {
    return zml.linalg.blas.rot(f32, n, x, incx, y, incy, c, s, .{}) catch unreachable;
}
export fn cblas_drot(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, c: f64, s: f64) void {
    return zml.linalg.blas.rot(f64, n, x, incx, y, incy, c, s, .{}) catch unreachable;
}
export fn cblas_csrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f32, s: f32) void {
    return zml.linalg.blas.rot(zml.cf32, n, x, incx, y, incy, c, s, .{}) catch unreachable;
}
export fn cblas_zdrot(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int, c: f64, s: f64) void {
    return zml.linalg.blas.rot(zml.cf64, n, x, incx, y, incy, c, s, .{}) catch unreachable;
}

export fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return zml.linalg.blas.rotg(f32, a, b, c, s, .{}) catch unreachable;
}
export fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return zml.linalg.blas.rotg(f64, a, b, c, s, .{}) catch unreachable;
}
export fn cblas_crotg(a: *anyopaque, b: *anyopaque, c: *f32, s: *anyopaque) void {
    return zml.linalg.blas.rotg(zml.cf32, a, b, c, s, .{}) catch unreachable;
}
export fn cblas_zrotg(a: *anyopaque, b: *anyopaque, c: *f64, s: *anyopaque) void {
    return zml.linalg.blas.rotg(zml.cf64, a, b, c, s, .{}) catch unreachable;
}

export fn cblas_srotm(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, param: [*c]f32) void {
    return zml.linalg.blas.rotm(f32, n, x, incx, y, incy, param, .{}) catch unreachable;
}
export fn cblas_drotm(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, param: [*c]f64) void {
    return zml.linalg.blas.rotm(f64, n, x, incx, y, incy, param, .{}) catch unreachable;
}

export fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void {
    return zml.linalg.blas.rotmg(f32, d1, d2, x1, y1, param, .{}) catch unreachable;
}
export fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void {
    return zml.linalg.blas.rotmg(f64, d1, d2, x1, y1, param, .{}) catch unreachable;
}

export fn cblas_sscal(n: c_int, alpha: f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.scal(f32, n, alpha, x, incx, .{}) catch unreachable;
}
export fn cblas_dscal(n: c_int, alpha: f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.scal(f64, n, alpha, x, incx, .{}) catch unreachable;
}
export fn cblas_cscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.scal(zml.cf32, n, alpha.*, x, incx, .{}) catch unreachable;
}
export fn cblas_zscal(n: c_int, alpha: *const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.scal(zml.cf64, n, alpha.*, x, incx, .{}) catch unreachable;
}
export fn cblas_csscal(n: c_int, alpha: f32, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.csscal(n, alpha, x, incx, .{}) catch unreachable;
}
export fn cblas_zdscal(n: c_int, alpha: f64, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.zdscal(n, alpha, x, incx, .{}) catch unreachable;
}

export fn cblas_sswap(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.swap(f32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_dswap(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.swap(f64, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_cswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.swap(zml.cf32, n, x, incx, y, incy, .{}) catch unreachable;
}
export fn cblas_zswap(n: c_int, x: [*c]anyopaque, incx: c_int, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.swap(zml.cf64, n, x, incx, y, incy, .{}) catch unreachable;
}

export fn cblas_isamax(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(f32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_idamax(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(f64, n, x, incx, .{}) catch unreachable;
}
export fn cblas_icamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(zml.cf32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_izamax(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamax(zml.cf64, n, x, incx, .{}) catch unreachable;
}

export fn cblas_isamin(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(f32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_idamin(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(f64, n, x, incx, .{}) catch unreachable;
}
export fn cblas_icamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(zml.cf32, n, x, incx, .{}) catch unreachable;
}
export fn cblas_izamin(n: c_int, x: [*c]const anyopaque, incx: c_int) c_uint {
    return zml.linalg.blas.iamin(zml.cf64, n, x, incx, .{}) catch unreachable;
}

export fn cblas_sgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.gbmv(f32, order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.gbmv(f64, order, transa, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_cgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gbmv(zml.cf32, order, transa, m, n, kl, ku, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gbmv(zml.cf64, order, transa, m, n, kl, ku, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_sgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.gemv(f32, order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.gemv(f64, order, transa, m, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_cgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gemv(zml.cf32, order, transa, m, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.gemv(zml.cf64, order, transa, m, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_sger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.ger(f32, order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_dger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.ger(f64, order, m, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_cgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.geru(zml.cf32, order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_zgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.geru(zml.cf64, order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_cgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.gerc(zml.cf32, order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_zgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.gerc(zml.cf64, order, m, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_chbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hbmv(zml.cf32, order, uplo, n, k, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zhbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hbmv(zml.cf64, order, uplo, n, k, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_chemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hemv(zml.cf32, order, uplo, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zhemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hemv(zml.cf64, order, uplo, n, alpha.*, a, lda, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_cher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her(zml.cf32, order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}
export fn cblas_zher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her(zml.cf64, order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}

export fn cblas_cher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her2(zml.cf32, order, uplo, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_zher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, a: [*c]anyopaque, lda: c_int) void {
    return zml.linalg.blas.her2(zml.cf64, order, uplo, n, alpha.*, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_chpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hpmv(zml.cf32, order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy, .{}) catch unreachable;
}
export fn cblas_zhpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, Ap: [*c]const anyopaque, x: [*c]const anyopaque, incx: c_int, beta: *const anyopaque, y: [*c]anyopaque, incy: c_int) void {
    return zml.linalg.blas.hpmv(zml.cf64, order, uplo, n, alpha.*, Ap, x, incx, beta.*, y, incy, .{}) catch unreachable;
}

export fn cblas_chpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr(zml.cf32, order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}
export fn cblas_zhpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const anyopaque, incx: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr(zml.cf64, order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}

export fn cblas_chpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr2(zml.cf32, order, uplo, n, alpha.*, x, incx, y, incy, Ap, .{}) catch unreachable;
}
export fn cblas_zhpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: [*c]const anyopaque, incx: c_int, y: [*c]const anyopaque, incy: c_int, Ap: [*c]anyopaque) void {
    return zml.linalg.blas.hpr2(zml.cf64, order, uplo, n, alpha.*, x, incx, y, incy, Ap, .{}) catch unreachable;
}

export fn cblas_ssbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.sbmv(f32, order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dsbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.sbmv(f64, order, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}

export fn cblas_sspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, Ap: [*c]const f32, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.spmv(f32, order, uplo, n, alpha, Ap, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, Ap: [*c]const f64, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.spmv(f64, order, uplo, n, alpha, Ap, x, incx, beta, y, incy, .{}) catch unreachable;
}

export fn cblas_sspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, Ap: [*c]f32) void {
    return zml.linalg.blas.spr(f32, order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}
export fn cblas_dspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, Ap: [*c]f64) void {
    return zml.linalg.blas.spr(f64, order, uplo, n, alpha, x, incx, Ap, .{}) catch unreachable;
}

export fn cblas_sspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, Ap: [*c]f32) void {
    return zml.linalg.blas.spr2(f32, order, uplo, n, alpha, x, incx, y, incy, Ap, .{}) catch unreachable;
}
export fn cblas_dspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, Ap: [*c]f64) void {
    return zml.linalg.blas.spr2(f64, order, uplo, n, alpha, x, incx, y, incy, Ap, .{}) catch unreachable;
}

export fn cblas_ssymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.symv(f32, order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}
export fn cblas_dsymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.symv(f64, order, uplo, n, alpha, a, lda, x, incx, beta, y, incy, .{}) catch unreachable;
}

export fn cblas_ssyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.syr(f32, order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}
export fn cblas_dsyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.syr(f64, order, uplo, n, alpha, x, incx, a, lda, .{}) catch unreachable;
}

export fn cblas_ssyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.syr2(f32, order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}
export fn cblas_dsyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.syr2(f64, order, uplo, n, alpha, x, incx, y, incy, a, lda, .{}) catch unreachable;
}

export fn cblas_stbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tbmv(f32, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tbmv(f64, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbmv(zml.cf32, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbmv(zml.cf64, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_stbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tbsv(f32, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tbsv(f64, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbsv(zml.cf32, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tbsv(zml.cf64, order, uplo, transa, diag, n, k, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_stpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tpmv(f32, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_dtpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tpmv(f64, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ctpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpmv(zml.cf32, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ztpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpmv(zml.cf64, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}

export fn cblas_stpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.tpsv(f32, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_dtpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.tpsv(f64, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ctpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpsv(zml.cf32, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}
export fn cblas_ztpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, Ap: [*c]const anyopaque, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.tpsv(zml.cf64, order, uplo, transa, diag, n, Ap, x, incx, .{}) catch unreachable;
}

export fn cblas_strmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.trmv(f32, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.trmv(f64, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trmv(zml.cf32, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trmv(zml.cf64, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_strsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.trsv(f32, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_dtrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.trsv(f64, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ctrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trsv(zml.cf32, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}
export fn cblas_ztrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const anyopaque, lda: c_int, x: [*c]anyopaque, incx: c_int) void {
    return zml.linalg.blas.trsv(zml.cf64, order, uplo, transa, diag, n, a, lda, x, incx, .{}) catch unreachable;
}

export fn cblas_sgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.gemm(f32, order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_dgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.gemm(f64, order, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_cgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.gemm(zml.cf32, order, transa, transb, m, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.gemm(zml.cf64, order, transa, transb, m, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_chemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.hemm(zml.cf32, order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zhemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.hemm(zml.cf64, order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_cherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: [*c]const anyopaque, lda: c_int, beta: f32, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.herk(zml.cf32, order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_zherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: [*c]const anyopaque, lda: c_int, beta: f64, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.herk(zml.cf64, order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}

export fn cblas_cher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.her2k(zml.cf32, order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.her2k(zml.cf64, order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_ssymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.symm(f32, order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_dsymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.symm(f64, order, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_csymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.symm(zml.cf32, order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zsymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.symm(zml.cf64, order, side, uplo, m, n, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_ssyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.syrk(f32, order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_dsyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.syrk(f64, order, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, .{}) catch unreachable;
}
export fn cblas_csyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syrk(zml.cf32, order, uplo, trans, n, k, alpha.*, a, lda, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zsyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syrk(zml.cf64, order, uplo, trans, n, k, alpha.*, a, lda, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_ssyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(f32, order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_dsyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(f64, order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_csyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(zml.cf32, order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}
export fn cblas_zsyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]const anyopaque, ldb: c_int, beta: *const anyopaque, c: [*c]anyopaque, ldc: c_int) void {
    return zml.linalg.blas.syr2k(zml.cf64, order, uplo, trans, n, k, alpha.*, a, lda, b, ldb, beta.*, c, ldc, .{}) catch unreachable;
}

export fn cblas_strmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]f32, ldb: c_int) void {
    return zml.linalg.blas.trmm(f32, order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_dtrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]f64, ldb: c_int) void {
    return zml.linalg.blas.trmm(f64, order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ctrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trmm(zml.cf32, order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ztrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trmm(zml.cf64, order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}

export fn cblas_strsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]f32, ldb: c_int) void {
    return zml.linalg.blas.trsm(f32, order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_dtrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]f64, ldb: c_int) void {
    return zml.linalg.blas.trsm(f64, order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ctrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trsm(zml.cf32, order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}
export fn cblas_ztrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: [*c]const anyopaque, lda: c_int, b: [*c]anyopaque, ldb: c_int) void {
    return zml.linalg.blas.trsm(zml.cf64, order, side, uplo, transa, diag, m, n, alpha.*, a, lda, b, ldb, .{}) catch unreachable;
}
