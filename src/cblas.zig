const std = @import("std");
const zsl = @import("zsl");

const cf32 = extern struct {
    re: f32,
    im: f32,
};

const cf64 = extern struct {
    re: f64,
    im: f64,
};

const CBLAS_LAYOUT = enum(c_int) {
    CblasRowMajor = 101,
    CblasColMajor = 102,

    fn to_zsl(self: CBLAS_LAYOUT) zsl.Layout {
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

    fn to_zsl(self: CBLAS_TRANSPOSE) zsl.linalg.Transpose {
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

    fn to_zsl(self: CBLAS_UPLO) zsl.Uplo {
        return switch (self) {
            .CblasUpper => .upper,
            .CblasLower => .lower,
        };
    }
};

const CBLAS_DIAG = enum(c_int) {
    CblasNonUnit = 131,
    CblasUnit = 132,

    fn to_zsl(self: CBLAS_DIAG) zsl.Diag {
        return switch (self) {
            .CblasNonUnit => .non_unit,
            .CblasUnit => .unit,
        };
    }
};

const CBLAS_SIDE = enum(c_int) {
    CblasLeft = 141,
    CblasRight = 142,

    fn to_zsl(self: CBLAS_SIDE) zsl.linalg.Side {
        return switch (self) {
            .CblasLeft => .left,
            .CblasRight => .right,
        };
    }
};

// Level 1
export fn cblas_sasum(n: isize, x: [*c]const f32, incx: isize) f32 {
    return zsl.linalg.blas.asum(zsl.numeric.cast(usize, n), @as([*]const f32, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}
export fn cblas_dasum(n: isize, x: [*c]const f64, incx: isize) f64 {
    return zsl.linalg.blas.asum(zsl.numeric.cast(usize, n), @as([*]const f64, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}
export fn cblas_scasum(n: isize, x: *const anyopaque, incx: isize) f32 {
    return zsl.linalg.blas.asum(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}
export fn cblas_dzasum(n: isize, x: *const anyopaque, incx: isize) f64 {
    return zsl.linalg.blas.asum(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}

export fn cblas_saxpy(n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]f32, incy: isize) void {
    return zsl.linalg.blas.axpy(zsl.numeric.cast(usize, n), alpha, @as([*]const f32, @ptrCast(@alignCast(x))), incx, @as([*]f32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_daxpy(n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]f64, incy: isize) void {
    return zsl.linalg.blas.axpy(zsl.numeric.cast(usize, n), alpha, @as([*]const f64, @ptrCast(@alignCast(x))), incx, @as([*]f64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_caxpy(n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *anyopaque, incy: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.axpy(zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_zaxpy(n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *anyopaque, incy: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.axpy(zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}

export fn cblas_scopy(n: isize, x: [*c]const f32, incx: isize, y: [*c]f32, incy: isize) void {
    return zsl.linalg.blas.copy(zsl.numeric.cast(usize, n), @as([*]const f32, @ptrCast(@alignCast(x))), incx, @as([*]f32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_dcopy(n: isize, x: [*c]const f64, incx: isize, y: [*c]f64, incy: isize) void {
    return zsl.linalg.blas.copy(zsl.numeric.cast(usize, n), @as([*]const f64, @ptrCast(@alignCast(x))), incx, @as([*]f64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_ccopy(n: isize, x: *const anyopaque, incx: isize, y: *anyopaque, incy: isize) void {
    return zsl.linalg.blas.copy(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_zcopy(n: isize, x: *const anyopaque, incx: isize, y: *anyopaque, incy: isize) void {
    return zsl.linalg.blas.copy(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}

export fn cblas_sdot(n: isize, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize) f32 {
    return zsl.linalg.blas.dot(zsl.numeric.cast(usize, n), @as([*]const f32, @ptrCast(@alignCast(x))), incx, @as([*]const f32, @ptrCast(@alignCast(y))), incy, .{}) catch 0.0;
}
export fn cblas_ddot(n: isize, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize) f64 {
    return zsl.linalg.blas.dot(zsl.numeric.cast(usize, n), @as([*]const f64, @ptrCast(@alignCast(x))), incx, @as([*]const f64, @ptrCast(@alignCast(y))), incy, .{}) catch 0.0;
}

export fn cblas_cdotc(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize) cf32 {
    const result: zsl.cf32 = zsl.linalg.blas.dotc(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch .{ .re = 0.0, .im = 0.0 };
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_zdotc(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize) cf64 {
    const result: zsl.cf64 = zsl.linalg.blas.dotc(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch .{ .re = 0.0, .im = 0.0 };
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_cdotc_sub(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, ret: *anyopaque) void {
    const ret_: *zsl.cf32 = @ptrCast(@alignCast(ret));
    return zsl.numeric.set(ret_, zsl.linalg.blas.dotc(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch zsl.cf32{ .re = 0.0, .im = 0.0 });
}
export fn cblas_zdotc_sub(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, ret: *anyopaque) void {
    const ret_: *zsl.cf64 = @ptrCast(@alignCast(ret));
    return zsl.numeric.set(ret_, zsl.linalg.blas.dotc(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch zsl.cf64{ .re = 0.0, .im = 0.0 });
}

export fn cblas_cdotu(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize) cf32 {
    const result: zsl.cf32 = zsl.linalg.blas.dot(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch .{ .re = 0.0, .im = 0.0 };
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_zdotu(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize) cf64 {
    const result: zsl.cf64 = zsl.linalg.blas.dot(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch .{ .re = 0.0, .im = 0.0 };
    return .{ .re = result.re, .im = result.im };
}
export fn cblas_cdotu_sub(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, ret: *anyopaque) void {
    const ret_: *zsl.cf32 = @ptrCast(@alignCast(ret));
    return zsl.numeric.set(ret_, zsl.linalg.blas.dot(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch zsl.cf32{ .re = 0.0, .im = 0.0 });
}
export fn cblas_zdotu_sub(n: isize, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, ret: *anyopaque) void {
    const ret_: *zsl.cf64 = @ptrCast(@alignCast(ret));
    return zsl.numeric.set(ret_, zsl.linalg.blas.dot(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch zsl.cf64{ .re = 0.0, .im = 0.0 });
}

export fn cblas_snrm2(n: isize, x: [*c]const f32, incx: isize) f32 {
    return zsl.linalg.blas.nrm2(zsl.numeric.cast(usize, n), @as([*]const f32, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}
export fn cblas_dnrm2(n: isize, x: [*c]const f64, incx: isize) f64 {
    return zsl.linalg.blas.nrm2(zsl.numeric.cast(usize, n), @as([*]const f64, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}
export fn cblas_scnrm2(n: isize, x: *const anyopaque, incx: isize) f32 {
    @setEvalBranchQuota(10_000);
    return zsl.linalg.blas.nrm2(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}
export fn cblas_dznrm2(n: isize, x: *const anyopaque, incx: isize) f64 {
    return zsl.linalg.blas.nrm2(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, .{}) catch 0.0;
}

export fn cblas_srot(n: isize, x: [*c]f32, incx: isize, y: [*c]f32, incy: isize, c: f32, s: f32) void {
    return zsl.linalg.blas.rot(zsl.numeric.cast(usize, n), @as([*]f32, @ptrCast(@alignCast(x))), incx, @as([*]f32, @ptrCast(@alignCast(y))), incy, c, s, .{}) catch {};
}
export fn cblas_drot(n: isize, x: [*c]f64, incx: isize, y: [*c]f64, incy: isize, c: f64, s: f64) void {
    return zsl.linalg.blas.rot(zsl.numeric.cast(usize, n), @as([*]f64, @ptrCast(@alignCast(x))), incx, @as([*]f64, @ptrCast(@alignCast(y))), incy, c, s, .{}) catch {};
}
export fn cblas_csrot(n: isize, x: *anyopaque, incx: isize, y: *anyopaque, incy: isize, c: f32, s: f32) void {
    return zsl.linalg.blas.rot(zsl.numeric.cast(usize, n), @as([*]zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf32, @ptrCast(@alignCast(y))), incy, c, s, .{}) catch {};
}
export fn cblas_zdrot(n: isize, x: *anyopaque, incx: isize, y: *anyopaque, incy: isize, c: f64, s: f64) void {
    return zsl.linalg.blas.rot(zsl.numeric.cast(usize, n), @as([*]zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf64, @ptrCast(@alignCast(y))), incy, c, s, .{}) catch {};
}

// export fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
//     return zsl.linalg.blas.srotg(a, b, c, s);
// }
// export fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
//     return zsl.linalg.blas.drotg(a, b, c, s);
// }
// export fn cblas_crotg(a: *anyopaque, b: *anyopaque, c: *f32, s: *anyopaque) void {
//     return zsl.linalg.blas.crotg(@ptrCast(@alignCast(a)), @ptrCast(@alignCast(b)), c, @ptrCast(@alignCast(s)));
// }
// export fn cblas_zrotg(a: *anyopaque, b: *anyopaque, c: *f64, s: *anyopaque) void {
//     return zsl.linalg.blas.zrotg(@ptrCast(@alignCast(a)), @ptrCast(@alignCast(b)), c, @ptrCast(@alignCast(s)));
// }

// export fn cblas_srotm(n: isize, x: [*c]f32, incx: isize, y: [*c]f32, incy: isize, param: [*c]f32) void {
//     return zsl.linalg.blas.srotm(n, x, incx, y, incy, param);
// }
// export fn cblas_drotm(n: isize, x: [*c]f64, incx: isize, y: [*c]f64, incy: isize, param: [*c]f64) void {
//     return zsl.linalg.blas.drotm(n, x, incx, y, incy, param);
// }

// export fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void {
//     return zsl.linalg.blas.srotmg(d1, d2, x1, y1, param);
// }
// export fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void {
//     return zsl.linalg.blas.drotmg(d1, d2, x1, y1, param);
// }

export fn cblas_sscal(n: isize, alpha: f32, x: [*c]f32, incx: isize) void {
    return zsl.linalg.blas.scal(zsl.numeric.cast(usize, n), alpha, @as([*]f32, @ptrCast(@alignCast(x))), incx, .{}) catch {};
}
export fn cblas_dscal(n: isize, alpha: f64, x: [*c]f64, incx: isize) void {
    return zsl.linalg.blas.scal(zsl.numeric.cast(usize, n), alpha, @as([*]f64, @ptrCast(@alignCast(x))), incx, .{}) catch {};
}
export fn cblas_cscal(n: isize, alpha: *const anyopaque, x: *anyopaque, incx: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.scal(zsl.numeric.cast(usize, n), alpha_.*, @as([*]zsl.cf32, @ptrCast(@alignCast(x))), incx, .{}) catch {};
}
export fn cblas_zscal(n: isize, alpha: *const anyopaque, x: *anyopaque, incx: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.scal(zsl.numeric.cast(usize, n), alpha_.*, @as([*]zsl.cf64, @ptrCast(@alignCast(x))), incx, .{}) catch {};
}
export fn cblas_csscal(n: isize, alpha: f32, x: *anyopaque, incx: isize) void {
    return zsl.linalg.blas.scal(zsl.numeric.cast(usize, n), alpha, @as([*]zsl.cf32, @ptrCast(@alignCast(x))), incx, .{}) catch {};
}
export fn cblas_zdscal(n: isize, alpha: f64, x: *anyopaque, incx: isize) void {
    return zsl.linalg.blas.scal(zsl.numeric.cast(usize, n), alpha, @as([*]zsl.cf64, @ptrCast(@alignCast(x))), incx, .{}) catch {};
}

export fn cblas_sswap(n: isize, x: [*c]f32, incx: isize, y: [*c]f32, incy: isize) void {
    return zsl.linalg.blas.swap(zsl.numeric.cast(usize, n), @as([*]f32, @ptrCast(@alignCast(x))), incx, @as([*]f32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_dswap(n: isize, x: [*c]f64, incx: isize, y: [*c]f64, incy: isize) void {
    return zsl.linalg.blas.swap(zsl.numeric.cast(usize, n), @as([*]f64, @ptrCast(@alignCast(x))), incx, @as([*]f64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_cswap(n: isize, x: *anyopaque, incx: isize, y: *anyopaque, incy: isize) void {
    return zsl.linalg.blas.swap(zsl.numeric.cast(usize, n), @as([*]zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_zswap(n: isize, x: *anyopaque, incx: isize, y: *anyopaque, incy: isize) void {
    return zsl.linalg.blas.swap(zsl.numeric.cast(usize, n), @as([*]zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}

export fn cblas_isamax(n: isize, x: [*c]const f32, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamax(zsl.numeric.cast(usize, n), @as([*]const f32, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}
export fn cblas_idamax(n: isize, x: [*c]const f64, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamax(zsl.numeric.cast(usize, n), @as([*]const f64, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}
export fn cblas_icamax(n: isize, x: *const anyopaque, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamax(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}
export fn cblas_izamax(n: isize, x: *const anyopaque, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamax(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}

export fn cblas_isamin(n: isize, x: [*c]const f32, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamin(zsl.numeric.cast(usize, n), @as([*]const f32, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}
export fn cblas_idamin(n: isize, x: [*c]const f64, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamin(zsl.numeric.cast(usize, n), @as([*]const f64, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}
export fn cblas_icamin(n: isize, x: *const anyopaque, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamin(zsl.numeric.cast(usize, n), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}
export fn cblas_izamin(n: isize, x: *const anyopaque, incx: isize) c_uint {
    return zsl.numeric.cast(c_uint, zsl.linalg.blas.iamin(zsl.numeric.cast(usize, n), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, .{}) catch 0);
}

// Level 2
export fn cblas_sgemv(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void {
    return zsl.linalg.blas.gemv(layout.to_zsl(), transa.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha, @as([*]const f32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const f32, @ptrCast(@alignCast(x))), zsl.numeric.cast(isize, incx), beta, @as([*]f32, @ptrCast(@alignCast(y))), zsl.numeric.cast(isize, incy), .{}) catch {};
}
export fn cblas_dgemv(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void {
    return zsl.linalg.blas.gemv(layout.to_zsl(), transa.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha, @as([*]const f64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const f64, @ptrCast(@alignCast(x))), zsl.numeric.cast(isize, incx), beta, @as([*]f64, @ptrCast(@alignCast(y))), zsl.numeric.cast(isize, incy), .{}) catch {};
}
export fn cblas_cgemv(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, x: *const anyopaque, incx: isize, beta: *const anyopaque, y: *anyopaque, incy: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
    return zsl.linalg.blas.gemv(layout.to_zsl(), transa.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), zsl.numeric.cast(isize, incx), beta_.*, @as([*]zsl.cf32, @ptrCast(@alignCast(y))), zsl.numeric.cast(isize, incy), .{}) catch {};
}
export fn cblas_zgemv(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, x: *const anyopaque, incx: isize, beta: *const anyopaque, y: *anyopaque, incy: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
    return zsl.linalg.blas.gemv(layout.to_zsl(), transa.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), zsl.numeric.cast(isize, incx), beta_.*, @as([*]zsl.cf64, @ptrCast(@alignCast(y))), zsl.numeric.cast(isize, incy), .{}) catch {};
}

export fn cblas_sger(layout: CBLAS_LAYOUT, m: isize, n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize, a: [*c]f32, lda: isize) void {
    return zsl.linalg.blas.ger(layout.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha, @as([*]const f32, @ptrCast(@alignCast(x))), incx, @as([*]const f32, @ptrCast(@alignCast(y))), incy, @as([*]f32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_dger(layout: CBLAS_LAYOUT, m: isize, n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize, a: [*c]f64, lda: isize) void {
    return zsl.linalg.blas.ger(layout.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha, @as([*]const f64, @ptrCast(@alignCast(x))), incx, @as([*]const f64, @ptrCast(@alignCast(y))), incy, @as([*]f64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

export fn cblas_cgeru(layout: CBLAS_LAYOUT, m: isize, n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, a: *anyopaque, lda: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.ger(layout.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, @as([*]zsl.cf32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_zgeru(layout: CBLAS_LAYOUT, m: isize, n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, a: *anyopaque, lda: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.ger(layout.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, @as([*]zsl.cf64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

export fn cblas_cgerc(layout: CBLAS_LAYOUT, m: isize, n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, a: *anyopaque, lda: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.gerc(layout.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, @as([*]zsl.cf32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_zgerc(layout: CBLAS_LAYOUT, m: isize, n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, a: *anyopaque, lda: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.gerc(layout.to_zsl(), zsl.numeric.cast(usize, m), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, @as([*]zsl.cf64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

export fn cblas_chemv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, x: *const anyopaque, incx: isize, beta: *const anyopaque, y: *anyopaque, incy: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
    return zsl.linalg.blas.hemv(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, beta_.*, @as([*]zsl.cf32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_zhemv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, x: *const anyopaque, incx: isize, beta: *const anyopaque, y: *anyopaque, incy: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
    return zsl.linalg.blas.hemv(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, beta_.*, @as([*]zsl.cf64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}

export fn cblas_cher(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f32, x: *const anyopaque, incx: isize, a: *anyopaque, lda: isize) void {
    return zsl.linalg.blas.her(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_zher(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f64, x: *const anyopaque, incx: isize, a: *anyopaque, lda: isize) void {
    return zsl.linalg.blas.her(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]zsl.cf64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

export fn cblas_cher2(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, a: *anyopaque, lda: isize) void {
    const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.her2(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf32, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf32, @ptrCast(@alignCast(y))), incy, @as([*]zsl.cf32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_zher2(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: *const anyopaque, x: *const anyopaque, incx: isize, y: *const anyopaque, incy: isize, a: *anyopaque, lda: isize) void {
    const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
    return zsl.linalg.blas.her2(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha_.*, @as([*]const zsl.cf64, @ptrCast(@alignCast(x))), incx, @as([*]const zsl.cf64, @ptrCast(@alignCast(y))), incy, @as([*]zsl.cf64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

export fn cblas_ssymv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f32, a: [*c]const f32, lda: isize, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void {
    return zsl.linalg.blas.symv(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const f32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const f32, @ptrCast(@alignCast(x))), incx, beta, @as([*]f32, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}
export fn cblas_dsymv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f64, a: [*c]const f64, lda: isize, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void {
    return zsl.linalg.blas.symv(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const f64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), @as([*]const f64, @ptrCast(@alignCast(x))), incx, beta, @as([*]f64, @ptrCast(@alignCast(y))), incy, .{}) catch {};
}

export fn cblas_ssyr(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f32, x: [*c]const f32, incx: isize, a: [*c]f32, lda: isize) void {
    return zsl.linalg.blas.syr(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const f32, @ptrCast(@alignCast(x))), incx, @as([*]f32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_dsyr(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f64, x: [*c]const f64, incx: isize, a: [*c]f64, lda: isize) void {
    return zsl.linalg.blas.syr(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const f64, @ptrCast(@alignCast(x))), incx, @as([*]f64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

export fn cblas_ssyr2(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize, a: [*c]f32, lda: isize) void {
    return zsl.linalg.blas.syr2(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const f32, @ptrCast(@alignCast(x))), incx, @as([*]const f32, @ptrCast(@alignCast(y))), incy, @as([*]f32, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}
export fn cblas_dsyr2(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize, a: [*c]f64, lda: isize) void {
    return zsl.linalg.blas.syr2(layout.to_zsl(), uplo.to_zsl(), zsl.numeric.cast(usize, n), alpha, @as([*]const f64, @ptrCast(@alignCast(x))), incx, @as([*]const f64, @ptrCast(@alignCast(y))), incy, @as([*]f64, @ptrCast(@alignCast(a))), zsl.numeric.cast(usize, lda), .{}) catch {};
}

// export fn cblas_strmv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: [*c]const f32, lda: isize, x: [*c]f32, incx: isize) void {
//     return zsl.linalg.blas.strmv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }
// export fn cblas_dtrmv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: [*c]const f64, lda: isize, x: [*c]f64, incx: isize) void {
//     return zsl.linalg.blas.dtrmv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }
// export fn cblas_ctrmv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: *const anyopaque, lda: isize, x: *anyopaque, incx: isize) void {
//     return zsl.linalg.blas.ctrmv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }
// export fn cblas_ztrmv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: *const anyopaque, lda: isize, x: *anyopaque, incx: isize) void {
//     return zsl.linalg.blas.ztrmv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }

// export fn cblas_strsv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: [*c]const f32, lda: isize, x: [*c]f32, incx: isize) void {
//     return zsl.linalg.blas.strsv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }
// export fn cblas_dtrsv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: [*c]const f64, lda: isize, x: [*c]f64, incx: isize) void {
//     return zsl.linalg.blas.dtrsv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }
// export fn cblas_ctrsv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: *const anyopaque, lda: isize, x: *anyopaque, incx: isize) void {
//     return zsl.linalg.blas.ctrsv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }
// export fn cblas_ztrsv(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, n: isize, a: *const anyopaque, lda: isize, x: *anyopaque, incx: isize) void {
//     return zsl.linalg.blas.ztrsv(layout.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), n, @ptrCast(@alignCast(a)), lda, x, incx);
// }

// // Level 3
// export fn cblas_sgemm(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: isize, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, beta: f32, c: [*c]f32, ldc: isize) void {
//     return zsl.linalg.blas.sgemm(layout.to_zsl(), transa.to_zsl(), transb.to_zsl(), m, n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_dgemm(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: isize, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, beta: f64, c: [*c]f64, ldc: isize) void {
//     return zsl.linalg.blas.dgemm(layout.to_zsl(), transa.to_zsl(), transb.to_zsl(), m, n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_cgemm(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: isize, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.cgemm(layout.to_zsl(), transa.to_zsl(), transb.to_zsl(), m, n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zgemm(layout: CBLAS_LAYOUT, transa: CBLAS_TRANSPOSE, transb: CBLAS_TRANSPOSE, m: isize, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.zgemm(layout.to_zsl(), transa.to_zsl(), transb.to_zsl(), m, n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_chemm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.chemm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zhemm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.zhemm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_cherk(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: f32, a: *const anyopaque, lda: isize, beta: f32, c: *anyopaque, ldc: isize) void {
//     return zsl.linalg.blas.cherk(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zherk(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: f64, a: *const anyopaque, lda: isize, beta: f64, c: *anyopaque, ldc: isize) void {
//     return zsl.linalg.blas.zherk(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_cher2k(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: f32, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     return zsl.linalg.blas.cher2k(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zher2k(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: f64, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     return zsl.linalg.blas.zher2k(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_ssymm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, beta: f32, c: [*c]f32, ldc: isize) void {
//     return zsl.linalg.blas.ssymm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_dsymm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, beta: f64, c: [*c]f64, ldc: isize) void {
//     return zsl.linalg.blas.dsymm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_csymm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.csymm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zsymm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.zsymm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_ssyrk(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, beta: f32, c: [*c]f32, ldc: isize) void {
//     return zsl.linalg.blas.ssyrk(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_dsyrk(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, beta: f64, c: [*c]f64, ldc: isize) void {
//     return zsl.linalg.blas.dsyrk(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha, @ptrCast(@alignCast(a)), lda, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_csyrk(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.csyrk(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zsyrk(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.zsyrk(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_ssyr2k(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, beta: f32, c: [*c]f32, ldc: isize) void {
//     return zsl.linalg.blas.ssyr2k(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_dsyr2k(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, beta: f64, c: [*c]f64, ldc: isize) void {
//     return zsl.linalg.blas.dsyr2k(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_csyr2k(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf32 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.csyr2k(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }
// export fn cblas_zsyr2k(layout: CBLAS_LAYOUT, uplo: CBLAS_UPLO, trans: CBLAS_TRANSPOSE, n: isize, k: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *const anyopaque, ldb: isize, beta: *const anyopaque, c: *anyopaque, ldc: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     const beta_: *const zsl.cf64 = @ptrCast(@alignCast(beta));
//     return zsl.linalg.blas.zsyr2k(layout.to_zsl(), uplo.to_zsl(), trans.to_zsl(), n, k, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb, beta_.*, @ptrCast(@alignCast(c)), ldc);
// }

// export fn cblas_strmm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) void {
//     return zsl.linalg.blas.strmm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
// export fn cblas_dtrmm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) void {
//     return zsl.linalg.blas.dtrmm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
// export fn cblas_ctrmm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *anyopaque, ldb: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     return zsl.linalg.blas.ctrmm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
// export fn cblas_ztrmm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *anyopaque, ldb: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     return zsl.linalg.blas.ztrmm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }

// export fn cblas_strsm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) void {
//     return zsl.linalg.blas.strsm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
// export fn cblas_dtrsm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) void {
//     return zsl.linalg.blas.dtrsm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
// export fn cblas_ctrsm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *anyopaque, ldb: isize) void {
//     const alpha_: *const zsl.cf32 = @ptrCast(@alignCast(alpha));
//     return zsl.linalg.blas.ctrsm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
// export fn cblas_ztrsm(layout: CBLAS_LAYOUT, side: CBLAS_SIDE, uplo: CBLAS_UPLO, transa: CBLAS_TRANSPOSE, diag: CBLAS_DIAG, m: isize, n: isize, alpha: *const anyopaque, a: *const anyopaque, lda: isize, b: *anyopaque, ldb: isize) void {
//     const alpha_: *const zsl.cf64 = @ptrCast(@alignCast(alpha));
//     return zsl.linalg.blas.ztrsm(layout.to_zsl(), side.to_zsl(), uplo.to_zsl(), transa.to_zsl(), diag.to_zsl(), m, n, alpha_.*, @ptrCast(@alignCast(a)), lda, @ptrCast(@alignCast(b)), ldb);
// }
