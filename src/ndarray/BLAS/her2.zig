const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Uplo = @import("../ndarray.zig").Uplo;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn her2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: scalar(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;
    const LENY = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.her does not support bool."),
        .BuiltinInt, .BuiltinFloat => @compileError("BLAS.her does not support int or float."),
        .Complex => {
            if (alpha.re == 0 and alpha.im == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.im * x[@intCast(jx)].re + alpha.re * x[@intCast(jx)].im);
                    const t1 = T.init(alpha.re * y[@intCast(jy)].re + alpha.im * y[@intCast(jy)].im, -alpha.im * y[@intCast(jy)].re + alpha.re * y[@intCast(jy)].im);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix = if (incx < 0) (-LENX + 1) * incx else 0;
                    var iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        A[@intCast(iaij)].re += t0.re * y[@intCast(iy)].re + t0.im * y[@intCast(iy)].im + t1.re * x[@intCast(ix)].re + t1.im * x[@intCast(ix)].im;
                        A[@intCast(iaij)].im += (t0.im * y[@intCast(iy)].re - t0.re * y[@intCast(iy)].im + t1.im * x[@intCast(ix)].re - t1.re * x[@intCast(ix)].im) * (-conj);

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    A[@intCast(iaij)].re += 2 * (t0.re * y[@intCast(jy)].re + t0.im * y[@intCast(jy)].im);
                    A[@intCast(iaij)].im = 0;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                    jy += incy;
                }
            } else {
                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.im * x[@intCast(jx)].re + alpha.re * x[@intCast(jx)].im);
                    const t1 = T.init(alpha.re * y[@intCast(jy)].re + alpha.im * y[@intCast(jy)].im, -alpha.im * y[@intCast(jy)].re + alpha.re * y[@intCast(jy)].im);
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    A[@intCast(jaj)].re += 2 * (t0.re * y[@intCast(jy)].re + t0.im * y[@intCast(jy)].im);
                    A[@intCast(jaj)].im = 0;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    var iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        A[@intCast(iaij)].re += t0.re * y[@intCast(iy)].re + t0.im * y[@intCast(iy)].im + t1.re * x[@intCast(ix)].re + t1.im * x[@intCast(ix)].im;
                        A[@intCast(iaij)].im += (t0.im * y[@intCast(iy)].re - t0.re * y[@intCast(iy)].im + t1.im * x[@intCast(ix)].re - t1.re * x[@intCast(ix)].im) * (-conj);

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.her only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "her2" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const alpha = Complex(f64).init(1, 1);

    const A = try a.alloc(Complex(f64), n * n);
    defer a.free(A);
    const x1 = try a.alloc(Complex(f64), n);
    defer a.free(x1);
    const y1 = try a.alloc(Complex(f64), n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 0),
        Complex(f64).init(2, -1),
        Complex(f64).init(3, 2),
        Complex(f64).init(4, -3),
        Complex(f64).init(5, 1),
        Complex(f64).init(2, 1),
        Complex(f64).init(6, 0),
        Complex(f64).init(7, -2),
        Complex(f64).init(8, 1),
        Complex(f64).init(9, -4),
        Complex(f64).init(3, -2),
        Complex(f64).init(7, 2),
        Complex(f64).init(10, 0),
        Complex(f64).init(11, -3),
        Complex(f64).init(12, 2),
        Complex(f64).init(4, 3),
        Complex(f64).init(8, -1),
        Complex(f64).init(11, 3),
        Complex(f64).init(13, 0),
        Complex(f64).init(14, -1),
        Complex(f64).init(5, -1),
        Complex(f64).init(9, 4),
        Complex(f64).init(12, -2),
        Complex(f64).init(14, 1),
        Complex(f64).init(15, 0),
    });
    @memcpy(x1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });
    @memcpy(y1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.her2(Complex(f64), .RowMajor, .Upper, n, alpha, x1.ptr, 1, y1.ptr, 1, A.ptr, n);

    try std.testing.expectEqual(11, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(21, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(29, A[2].re);
    try std.testing.expectEqual(-6, A[2].im);
    try std.testing.expectEqual(39, A[3].re);
    try std.testing.expectEqual(-8, A[3].im);
    try std.testing.expectEqual(47, A[4].re);
    try std.testing.expectEqual(-15, A[4].im);
    try std.testing.expectEqual(2, A[5].re);
    try std.testing.expectEqual(1, A[5].im);
    try std.testing.expectEqual(42, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(50, A[7].re);
    try std.testing.expectEqual(-29, A[7].im);
    try std.testing.expectEqual(68, A[8].re);
    try std.testing.expectEqual(-23, A[8].im);
    try std.testing.expectEqual(76, A[9].re);
    try std.testing.expectEqual(-55, A[9].im);
    try std.testing.expectEqual(3, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(7, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(52, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(70, A[13].re);
    try std.testing.expectEqual(16, A[13].im);
    try std.testing.expectEqual(70, A[14].re);
    try std.testing.expectEqual(-6, A[14].im);
    try std.testing.expectEqual(4, A[15].re);
    try std.testing.expectEqual(3, A[15].im);
    try std.testing.expectEqual(8, A[16].re);
    try std.testing.expectEqual(-1, A[16].im);
    try std.testing.expectEqual(11, A[17].re);
    try std.testing.expectEqual(3, A[17].im);
    try std.testing.expectEqual(97, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(97, A[19].re);
    try std.testing.expectEqual(-44, A[19].im);
    try std.testing.expectEqual(5, A[20].re);
    try std.testing.expectEqual(-1, A[20].im);
    try std.testing.expectEqual(9, A[21].re);
    try std.testing.expectEqual(4, A[21].im);
    try std.testing.expectEqual(12, A[22].re);
    try std.testing.expectEqual(-2, A[22].im);
    try std.testing.expectEqual(14, A[23].re);
    try std.testing.expectEqual(1, A[23].im);
    try std.testing.expectEqual(89, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);
    const y2 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(7, 8),
        Complex(f64).init(5, 6),
        Complex(f64).init(3, 4),
        Complex(f64).init(1, 2),
    });
    @memcpy(y2.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });

    BLAS.her2(Complex(f64), .ColumnMajor, .Upper, n, alpha, x2.ptr, -1, y2.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(21, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(21, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(29, A[2].re);
    try std.testing.expectEqual(-6, A[2].im);
    try std.testing.expectEqual(39, A[3].re);
    try std.testing.expectEqual(-8, A[3].im);
    try std.testing.expectEqual(47, A[4].re);
    try std.testing.expectEqual(-15, A[4].im);
    try std.testing.expectEqual(21, A[5].re);
    try std.testing.expectEqual(4, A[5].im);
    try std.testing.expectEqual(78, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(50, A[7].re);
    try std.testing.expectEqual(-29, A[7].im);
    try std.testing.expectEqual(68, A[8].re);
    try std.testing.expectEqual(-23, A[8].im);
    try std.testing.expectEqual(76, A[9].re);
    try std.testing.expectEqual(-55, A[9].im);
    try std.testing.expectEqual(29, A[10].re);
    try std.testing.expectEqual(-10, A[10].im);
    try std.testing.expectEqual(50, A[11].re);
    try std.testing.expectEqual(-25, A[11].im);
    try std.testing.expectEqual(94, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(70, A[13].re);
    try std.testing.expectEqual(16, A[13].im);
    try std.testing.expectEqual(70, A[14].re);
    try std.testing.expectEqual(-6, A[14].im);
    try std.testing.expectEqual(39, A[15].re);
    try std.testing.expectEqual(-2, A[15].im);
    try std.testing.expectEqual(68, A[16].re);
    try std.testing.expectEqual(-25, A[16].im);
    try std.testing.expectEqual(70, A[17].re);
    try std.testing.expectEqual(22, A[17].im);
    try std.testing.expectEqual(181, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(97, A[19].re);
    try std.testing.expectEqual(-44, A[19].im);
    try std.testing.expectEqual(47, A[20].re);
    try std.testing.expectEqual(-17, A[20].im);
    try std.testing.expectEqual(76, A[21].re);
    try std.testing.expectEqual(-47, A[21].im);
    try std.testing.expectEqual(70, A[22].re);
    try std.testing.expectEqual(-10, A[22].im);
    try std.testing.expectEqual(97, A[23].re);
    try std.testing.expectEqual(-42, A[23].im);
    try std.testing.expectEqual(163, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x3 = try a.alloc(Complex(f64), n);
    defer a.free(x3);
    const y3 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });
    @memcpy(y3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });

    BLAS.her2(Complex(f64), .RowMajor, .Lower, n, alpha, x3.ptr, 1, y3.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(31, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(21, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(29, A[2].re);
    try std.testing.expectEqual(-6, A[2].im);
    try std.testing.expectEqual(39, A[3].re);
    try std.testing.expectEqual(-8, A[3].im);
    try std.testing.expectEqual(47, A[4].re);
    try std.testing.expectEqual(-15, A[4].im);
    try std.testing.expectEqual(40, A[5].re);
    try std.testing.expectEqual(1, A[5].im);
    try std.testing.expectEqual(114, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(50, A[7].re);
    try std.testing.expectEqual(-29, A[7].im);
    try std.testing.expectEqual(68, A[8].re);
    try std.testing.expectEqual(-23, A[8].im);
    try std.testing.expectEqual(76, A[9].re);
    try std.testing.expectEqual(-55, A[9].im);
    try std.testing.expectEqual(55, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(93, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(136, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(70, A[13].re);
    try std.testing.expectEqual(16, A[13].im);
    try std.testing.expectEqual(70, A[14].re);
    try std.testing.expectEqual(-6, A[14].im);
    try std.testing.expectEqual(74, A[15].re);
    try std.testing.expectEqual(3, A[15].im);
    try std.testing.expectEqual(128, A[16].re);
    try std.testing.expectEqual(-1, A[16].im);
    try std.testing.expectEqual(129, A[17].re);
    try std.testing.expectEqual(3, A[17].im);
    try std.testing.expectEqual(265, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(97, A[19].re);
    try std.testing.expectEqual(-44, A[19].im);
    try std.testing.expectEqual(89, A[20].re);
    try std.testing.expectEqual(-1, A[20].im);
    try std.testing.expectEqual(143, A[21].re);
    try std.testing.expectEqual(4, A[21].im);
    try std.testing.expectEqual(128, A[22].re);
    try std.testing.expectEqual(-2, A[22].im);
    try std.testing.expectEqual(180, A[23].re);
    try std.testing.expectEqual(1, A[23].im);
    try std.testing.expectEqual(237, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x4 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(Complex(f64), n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y4.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(1, 2),
    });

    BLAS.her2(Complex(f64), .ColumnMajor, .Lower, n, alpha, x4.ptr, -2, y4.ptr, -1, A.ptr, n);

    try std.testing.expectEqual(41, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(40, A[1].re);
    try std.testing.expectEqual(-1, A[1].im);
    try std.testing.expectEqual(55, A[2].re);
    try std.testing.expectEqual(2, A[2].im);
    try std.testing.expectEqual(74, A[3].re);
    try std.testing.expectEqual(-3, A[3].im);
    try std.testing.expectEqual(89, A[4].re);
    try std.testing.expectEqual(1, A[4].im);
    try std.testing.expectEqual(40, A[5].re);
    try std.testing.expectEqual(1, A[5].im);
    try std.testing.expectEqual(150, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(93, A[7].re);
    try std.testing.expectEqual(-2, A[7].im);
    try std.testing.expectEqual(128, A[8].re);
    try std.testing.expectEqual(1, A[8].im);
    try std.testing.expectEqual(143, A[9].re);
    try std.testing.expectEqual(-4, A[9].im);
    try std.testing.expectEqual(55, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(93, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(178, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(129, A[13].re);
    try std.testing.expectEqual(-3, A[13].im);
    try std.testing.expectEqual(128, A[14].re);
    try std.testing.expectEqual(2, A[14].im);
    try std.testing.expectEqual(74, A[15].re);
    try std.testing.expectEqual(3, A[15].im);
    try std.testing.expectEqual(128, A[16].re);
    try std.testing.expectEqual(-1, A[16].im);
    try std.testing.expectEqual(129, A[17].re);
    try std.testing.expectEqual(3, A[17].im);
    try std.testing.expectEqual(349, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(180, A[19].re);
    try std.testing.expectEqual(-1, A[19].im);
    try std.testing.expectEqual(89, A[20].re);
    try std.testing.expectEqual(-1, A[20].im);
    try std.testing.expectEqual(143, A[21].re);
    try std.testing.expectEqual(4, A[21].im);
    try std.testing.expectEqual(128, A[22].re);
    try std.testing.expectEqual(-2, A[22].im);
    try std.testing.expectEqual(180, A[23].re);
    try std.testing.expectEqual(1, A[23].im);
    try std.testing.expectEqual(311, A[24].re);
    try std.testing.expectEqual(0, A[24].im);
}
