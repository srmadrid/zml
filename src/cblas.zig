const std = @import("std");
const zml = @import("zml");

const cf32 = extern struct {
    re: f32,
    im: f32,
};

const cf64 = extern struct {
    re: f64,
    im: f64,
};

const CBLAS_ORDER = enum(c_int) {
    CblasRowMajor = 101,
    CblasColMajor = 102,

    fn to_zml(self: CBLAS_ORDER) zml.linalg.Order {
        return switch (self) {
            .CblasRowMajor => .row_major,
            .CblasColMajor => .col_major,
        };
    }
};

const CBLAS_TRANSPOSE = enum(c_int) {
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113,
    CblasConjNoTrans = 114,

    fn to_zml(self: CBLAS_TRANSPOSE) zml.linalg.Transpose {
        return switch (self) {
            .CblasNoTrans => .no_trans,
            .CblasTrans => .trans,
            .CblasConjTrans => .conj_trans,
            .CblasConjNoTrans => .conj_no_trans,
        };
    }
};

const CBLAS_UPLO = enum(c_int) {
    CblasUpper = 121,
    CblasLower = 122,

    fn to_zml(self: CBLAS_UPLO) zml.linalg.Uplo {
        return switch (self) {
            .CblasUpper => .upper,
            .CblasLower => .lower,
        };
    }
};

const CBLAS_DIAG = enum(c_int) {
    CblasNonUnit = 131,
    CblasUnit = 132,

    fn to_zml(self: CBLAS_DIAG) zml.linalg.Diag {
        return switch (self) {
            .CblasNonUnit => .non_unit,
            .CblasUnit => .unit,
        };
    }
};

const CBLAS_SIDE = enum(c_int) {
    CblasLeft = 141,
    CblasRight = 142,

    fn to_zml(self: CBLAS_SIDE) zml.linalg.Side {
        return switch (self) {
            .CblasLeft => .left,
            .CblasRight => .right,
        };
    }
};

// Level 1
export fn cblas_sasum(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.linalg.blas.sasum(n, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dasum(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.linalg.blas.dasum(n, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_scasum(n: c_int, x: *const anyopaque, incx: c_int) f32 {
    return zml.linalg.blas.scasum(n, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dzasum(n: c_int, x: *const anyopaque, incx: c_int) f64 {
    return zml.linalg.blas.dzasum(n, @ptrCast(@alignCast(x)), incx);
}

export fn cblas_saxpy(n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.saxpy(n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_daxpy(n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.daxpy(n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_caxpy(n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.caxpy(n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zaxpy(n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zaxpy(n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_scopy(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.scopy(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dcopy(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dcopy(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_ccopy(n: c_int, x: *const anyopaque, incx: c_int, y: *anyopaque, incy: c_int) void {
    return zml.linalg.blas.ccopy(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zcopy(n: c_int, x: *const anyopaque, incx: c_int, y: *anyopaque, incy: c_int) void {
    return zml.linalg.blas.zcopy(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_sdot(n: c_int, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int) f32 {
    return zml.linalg.blas.sdot(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_ddot(n: c_int, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int) f64 {
    return zml.linalg.blas.ddot(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_cdotc(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int) cf32 {
    const result: zml.cf32 = zml.linalg.blas.cdotc(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_zdotc(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int) cf64 {
    const result: zml.cf64 = zml.linalg.blas.zdotc(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_cdotc_sub(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.cdotc_sub(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ret)));
}
export fn cblas_zdotc_sub(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.zdotc_sub(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ret)));
}

export fn cblas_cdotu(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int) cf32 {
    const result: zml.cf32 = zml.linalg.blas.cdotu(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_zdotu(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int) cf64 {
    const result: zml.cf64 = zml.linalg.blas.zdotu(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_cdotu_sub(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.cdotu_sub(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ret)));
}
export fn cblas_zdotu_sub(n: c_int, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, ret: *anyopaque) void {
    return zml.linalg.blas.zdotu_sub(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ret)));
}

export fn cblas_snrm2(n: c_int, x: [*c]const f32, incx: c_int) f32 {
    return zml.linalg.blas.snrm2(n, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dnrm2(n: c_int, x: [*c]const f64, incx: c_int) f64 {
    return zml.linalg.blas.dnrm2(n, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_scnrm2(n: c_int, x: *const anyopaque, incx: c_int) f32 {
    return zml.linalg.blas.scnrm2(n, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dznrm2(n: c_int, x: *const anyopaque, incx: c_int) f64 {
    return zml.linalg.blas.dznrm2(n, @ptrCast(@alignCast(x)), incx);
}

export fn cblas_srot(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, c: f32, s: f32) void {
    return zml.linalg.blas.srot(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, c, s);
}
export fn cblas_drot(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, c: f64, s: f64) void {
    return zml.linalg.blas.drot(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, c, s);
}
export fn cblas_csrot(n: c_int, x: *anyopaque, incx: c_int, y: *anyopaque, incy: c_int, c: f32, s: f32) void {
    return zml.linalg.blas.csrot(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, c, s);
}
export fn cblas_zdrot(n: c_int, x: *anyopaque, incx: c_int, y: *anyopaque, incy: c_int, c: f64, s: f64) void {
    return zml.linalg.blas.zdrot(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, c, s);
}

export fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return zml.linalg.blas.srotg(a, b, c, s);
}
export fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return zml.linalg.blas.drotg(a, b, c, s);
}
export fn cblas_crotg(a: *anyopaque, b: *anyopaque, c: *f32, s: *anyopaque) void {
    return zml.linalg.blas.crotg(@ptrCast(@alignCast(a)), @ptrCast(@alignCast(b)), c, @ptrCast(@alignCast(s)));
}
export fn cblas_zrotg(a: *anyopaque, b: *anyopaque, c: *f64, s: *anyopaque) void {
    return zml.linalg.blas.zrotg(@ptrCast(@alignCast(a)), @ptrCast(@alignCast(b)), c, @ptrCast(@alignCast(s)));
}

export fn cblas_srotm(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int, param: [*c]f32) void {
    return zml.linalg.blas.srotm(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, param);
}
export fn cblas_drotm(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int, param: [*c]f64) void {
    return zml.linalg.blas.drotm(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, param);
}

export fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void {
    return zml.linalg.blas.srotmg(d1, d2, x1, y1, param);
}
export fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void {
    return zml.linalg.blas.drotmg(d1, d2, x1, y1, param);
}

export fn cblas_sscal(n: c_int, alpha: f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.sscal(n, alpha, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dscal(n: c_int, alpha: f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dscal(n, alpha, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_cscal(n: c_int, alpha: *const anyopaque, x: *anyopaque, incx: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.cscal(n, alpha_.*, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_zscal(n: c_int, alpha: *const anyopaque, x: *anyopaque, incx: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zscal(n, alpha_.*, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_csscal(n: c_int, alpha: f32, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.csscal(n, alpha, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_zdscal(n: c_int, alpha: f64, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.zdscal(n, alpha, @ptrCast(@alignCast(x)), incx);
}

export fn cblas_sswap(n: c_int, x: [*c]f32, incx: c_int, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.sswap(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dswap(n: c_int, x: [*c]f64, incx: c_int, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dswap(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_cswap(n: c_int, x: *anyopaque, incx: c_int, y: *anyopaque, incy: c_int) void {
    return zml.linalg.blas.cswap(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zswap(n: c_int, x: *anyopaque, incx: c_int, y: *anyopaque, incy: c_int) void {
    return zml.linalg.blas.zswap(n, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_isamax(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.isamax(n, @ptrCast(@alignCast(x)), incx));
}
export fn cblas_idamax(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.idamax(n, @ptrCast(@alignCast(x)), incx));
}
export fn cblas_icamax(n: c_int, x: *const anyopaque, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.icamax(n, @ptrCast(@alignCast(x)), incx));
}
export fn cblas_izamax(n: c_int, x: *const anyopaque, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.izamax(n, @ptrCast(@alignCast(x)), incx));
}

export fn cblas_isamin(n: c_int, x: [*c]const f32, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.isamin(n, @ptrCast(@alignCast(x)), incx));
}
export fn cblas_idamin(n: c_int, x: [*c]const f64, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.idamin(n, @ptrCast(@alignCast(x)), incx));
}
export fn cblas_icamin(n: c_int, x: *const anyopaque, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.icamin(n, @ptrCast(@alignCast(x)), incx));
}
export fn cblas_izamin(n: c_int, x: *const anyopaque, incx: c_int) c_uint {
    return zml.scast(c_uint, zml.linalg.blas.izamin(n, @ptrCast(@alignCast(x)), incx));
}

// Level 2
export fn cblas_sgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.sgbmv(order.to_zml(), transa.to_zml(), m, n, kl, ku, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dgbmv(order.to_zml(), transa.to_zml(), m, n, kl, ku, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_cgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.cgbmv(order.to_zml(), transa.to_zml(), m, n, kl, ku, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zgbmv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, kl: c_int, ku: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zgbmv(order.to_zml(), transa.to_zml(), m, n, kl, ku, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_sgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.sgemv(order.to_zml(), transa.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dgemv(order.to_zml(), transa.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_cgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.cgemv(order.to_zml(), transa.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zgemv(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zgemv(order.to_zml(), transa.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_sger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.sger(order.to_zml(), m, n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_dger(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.dger(order.to_zml(), m, n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_cgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, a: *anyopaque, lda: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.cgeru(order.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_zgeru(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, a: *anyopaque, lda: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zgeru(order.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_cgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, a: *anyopaque, lda: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.cgerc(order.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_zgerc(order: CBLAS_ORDER, m: c_int, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, a: *anyopaque, lda: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zgerc(order.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_chbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.chbmv(order.to_zml(), uplo.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zhbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zhbmv(order.to_zml(), uplo.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_chemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.chemv(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zhemv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zhemv(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_cher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: *const anyopaque, incx: c_int, a: *anyopaque, lda: c_int) void {
    return zml.linalg.blas.cher(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_zher(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: *const anyopaque, incx: c_int, a: *anyopaque, lda: c_int) void {
    return zml.linalg.blas.zher(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_cher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, a: *anyopaque, lda: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.cher2(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_zher2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, a: *anyopaque, lda: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zher2(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_chpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, ap: *const anyopaque, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.chpmv(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_zhpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, ap: *const anyopaque, x: *const anyopaque, incx: c_int, beta: *const anyopaque, y: *anyopaque, incy: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zhpmv(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx, beta_.*, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_chpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: *const anyopaque, incx: c_int, ap: *anyopaque) void {
    return zml.linalg.blas.chpr(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(ap)));
}
export fn cblas_zhpr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: *const anyopaque, incx: c_int, ap: *anyopaque) void {
    return zml.linalg.blas.zhpr(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(ap)));
}

export fn cblas_chpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, ap: *anyopaque) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.chpr2(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ap)));
}
export fn cblas_zhpr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: *const anyopaque, x: *const anyopaque, incx: c_int, y: *const anyopaque, incy: c_int, ap: *anyopaque) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zhpr2(order.to_zml(), uplo.to_zml(), n, alpha_.*, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ap)));
}

export fn cblas_ssbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.ssbmv(order.to_zml(), uplo.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dsbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dsbmv(order.to_zml(), uplo.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_sspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, ap: [*c]const f32, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.sspmv(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dspmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, ap: [*c]const f64, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dspmv(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_sspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, ap: [*c]f32) void {
    return zml.linalg.blas.sspr(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(ap)));
}
export fn cblas_dspr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, ap: [*c]f64) void {
    return zml.linalg.blas.dspr(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(ap)));
}

export fn cblas_sspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, ap: [*c]f32) void {
    return zml.linalg.blas.sspr2(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ap)));
}
export fn cblas_dspr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, ap: [*c]f64) void {
    return zml.linalg.blas.dspr2(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(ap)));
}

export fn cblas_ssymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, x: [*c]const f32, incx: c_int, beta: f32, y: [*c]f32, incy: c_int) void {
    return zml.linalg.blas.ssymv(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}
export fn cblas_dsymv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, x: [*c]const f64, incx: c_int, beta: f64, y: [*c]f64, incy: c_int) void {
    return zml.linalg.blas.dsymv(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx, beta, @ptrCast(@alignCast(y)), incy);
}

export fn cblas_ssyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.ssyr(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_dsyr(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.dsyr(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_ssyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f32, x: [*c]const f32, incx: c_int, y: [*c]const f32, incy: c_int, a: [*c]f32, lda: c_int) void {
    return zml.linalg.blas.ssyr2(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}
export fn cblas_dsyr2(order: CBLAS_ORDER, uplo: CBLAS_UPLO, n: c_int, alpha: f64, x: [*c]const f64, incx: c_int, y: [*c]const f64, incy: c_int, a: [*c]f64, lda: c_int) void {
    return zml.linalg.blas.dsyr2(order.to_zml(), uplo.to_zml(), n, alpha, @ptrCast(@alignCast(x)), incx, @ptrCast(@alignCast(y)), incy, @ptrCast(@alignCast(a)), lda);
}

export fn cblas_stbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.stbmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dtbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dtbmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ctbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ctbmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ztbmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ztbmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}

export fn cblas_stbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.stbsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dtbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dtbsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ctbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ctbsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ztbsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, k: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ztbsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, k, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}

export fn cblas_stpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.stpmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dtpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dtpmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ctpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: *const anyopaque, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ctpmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ztpmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: *const anyopaque, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ztpmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}

export fn cblas_stpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: [*c]const f32, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.stpsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dtpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: [*c]const f64, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dtpsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ctpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: *const anyopaque, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ctpsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ztpsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, ap: *const anyopaque, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ztpsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(ap)), @ptrCast(@alignCast(x)), incx);
}

export fn cblas_strmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.strmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dtrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dtrmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ctrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ctrmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ztrmv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ztrmv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}

export fn cblas_strsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f32, lda: c_int, x: [*c]f32, incx: c_int) void {
    return zml.linalg.blas.strsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_dtrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: [*c]const f64, lda: c_int, x: [*c]f64, incx: c_int) void {
    return zml.linalg.blas.dtrsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ctrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ctrsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}
export fn cblas_ztrsv(order: CBLAS_ORDER, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: c_int, a: *const anyopaque, lda: c_int, x: *anyopaque, incx: c_int) void {
    return zml.linalg.blas.ztrsv(order.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), n, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(x)), incx);
}

// Level 3
export fn cblas_sgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.sgemm(order.to_zml(), transa.to_zml(), transb.to_zml(), m, n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_dgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.dgemm(order.to_zml(), transa.to_zml(), transb.to_zml(), m, n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_cgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.cgemm(order.to_zml(), transa.to_zml(), transb.to_zml(), m, n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zgemm(order: CBLAS_ORDER, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: c_int, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zgemm(order.to_zml(), transa.to_zml(), transb.to_zml(), m, n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_chemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.chemm(order.to_zml(), side.to_zml(), uplo.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zhemm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zhemm(order.to_zml(), side.to_zml(), uplo.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_cherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: *const anyopaque, lda: c_int, beta: f32, c: *anyopaque, ldc: c_int) void {
    return zml.linalg.blas.cherk(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zherk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: *const anyopaque, lda: c_int, beta: f64, c: *anyopaque, ldc: c_int) void {
    return zml.linalg.blas.zherk(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_cher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: f32, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.cher2k(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zher2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: f64, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.zher2k(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_ssymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.ssymm(order.to_zml(), side.to_zml(), uplo.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_dsymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.dsymm(order.to_zml(), side.to_zml(), uplo.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_csymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.csymm(order.to_zml(), side.to_zml(), uplo.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zsymm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zsymm(order.to_zml(), side.to_zml(), uplo.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_ssyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.ssyrk(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_dsyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.dsyrk(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_csyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.csyrk(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, beta_.*, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zsyrk(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zsyrk(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, beta_.*, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_ssyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]const f32, ldb: c_int, beta: f32, c: [*c]f32, ldc: c_int) void {
    return zml.linalg.blas.ssyr2k(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_dsyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]const f64, ldb: c_int, beta: f64, c: [*c]f64, ldc: c_int) void {
    return zml.linalg.blas.dsyr2k(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_csyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf32 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.csyr2k(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}
export fn cblas_zsyr2k(order: CBLAS_ORDER, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: c_int, k: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *const anyopaque, ldb: c_int, beta: *const anyopaque, c: *anyopaque, ldc: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zml.cf64 = @ptrCast(@alignCast(beta));
    return zml.linalg.blas.zsyr2k(order.to_zml(), uplo.to_zml(), trans.to_zml(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
}

export fn cblas_strmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]f32, ldb: c_int) void {
    return zml.linalg.blas.strmm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
export fn cblas_dtrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]f64, ldb: c_int) void {
    return zml.linalg.blas.dtrmm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
export fn cblas_ctrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *anyopaque, ldb: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.ctrmm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
export fn cblas_ztrmm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *anyopaque, ldb: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.ztrmm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}

export fn cblas_strsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f32, a: [*c]const f32, lda: c_int, b: [*c]f32, ldb: c_int) void {
    return zml.linalg.blas.strsm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
export fn cblas_dtrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: f64, a: [*c]const f64, lda: c_int, b: [*c]f64, ldb: c_int) void {
    return zml.linalg.blas.dtrsm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
export fn cblas_ctrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *anyopaque, ldb: c_int) void {
    const alpha_: *const zml.cf32 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.ctrsm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
export fn cblas_ztrsm(order: CBLAS_ORDER, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: c_int, n: c_int, alpha: *const anyopaque, a: *const anyopaque, lda: c_int, b: *anyopaque, ldb: c_int) void {
    const alpha_: *const zml.cf64 = @ptrCast(@alignCast(alpha));
    return zml.linalg.blas.ztrsm(order.to_zml(), side.to_zml(), uplo.to_zml(), transa.to_zml(), diag.to_zml(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
}
