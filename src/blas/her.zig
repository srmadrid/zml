const std = @import("std");
const core = @import("../core/core.zig");
const Order = @import("../ndarray/ndarray.zig").Order;
const Uplo = @import("../ndarray/ndarray.zig").Uplo;
const blas = @import("blas.zig");

const scalar = core.supported.scalar;

pub inline fn her(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: scalar(T), x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
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

    switch (supported) {
        .BuiltinBool => @compileError("blas.her does not support bool."),
        .BuiltinInt, .BuiltinFloat => @compileError("blas.her does not support int or float."),
        .Complex => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = T.init(alpha * x[@intCast(jx)].re, alpha * x[@intCast(jx)].im * (-conj));

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < j) {
                        A[@intCast(iaij)].re += t0.re * x[@intCast(ix)].re - t0.im * x[@intCast(ix)].im * conj;
                        A[@intCast(iaij)].im += t0.im * x[@intCast(ix)].re + t0.re * x[@intCast(ix)].im * conj;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    A[@intCast(iaij)].re += t0.re * x[@intCast(jx)].re - t0.im * x[@intCast(jx)].im * conj;
                    A[@intCast(iaij)].im = 0;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else {
                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = T.init(alpha * x[@intCast(jx)].re, alpha * x[@intCast(jx)].im * (-conj));
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    A[@intCast(jaj)].re += t0.re * x[@intCast(jx)].re - t0.im * x[@intCast(jx)].im * conj;
                    A[@intCast(jaj)].im = 0;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    while (i < I1) {
                        A[@intCast(iaij)].re += t0.re * x[@intCast(ix)].re - t0.im * x[@intCast(ix)].im * conj;
                        A[@intCast(iaij)].im += t0.im * x[@intCast(ix)].re + t0.re * x[@intCast(ix)].im * conj;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("blas.her only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "her" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const alpha: f64 = 2;

    const A = try a.alloc(Complex(f64), n * n);
    defer a.free(A);
    const x1 = try a.alloc(Complex(f64), n);
    defer a.free(x1);

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

    blas.her(Complex(f64), .RowMajor, .Upper, n, alpha, x1.ptr, 1, A.ptr, n);

    try std.testing.expectEqual(11, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(24, A[1].re);
    try std.testing.expectEqual(3, A[1].im);
    try std.testing.expectEqual(37, A[2].re);
    try std.testing.expectEqual(10, A[2].im);
    try std.testing.expectEqual(50, A[3].re);
    try std.testing.expectEqual(9, A[3].im);
    try std.testing.expectEqual(63, A[4].re);
    try std.testing.expectEqual(17, A[4].im);
    try std.testing.expectEqual(2, A[5].re);
    try std.testing.expectEqual(1, A[5].im);
    try std.testing.expectEqual(56, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(85, A[7].re);
    try std.testing.expectEqual(2, A[7].im);
    try std.testing.expectEqual(114, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(143, A[9].re);
    try std.testing.expectEqual(8, A[9].im);
    try std.testing.expectEqual(3, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(7, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(132, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(177, A[13].re);
    try std.testing.expectEqual(1, A[13].im);
    try std.testing.expectEqual(222, A[14].re);
    try std.testing.expectEqual(10, A[14].im);
    try std.testing.expectEqual(4, A[15].re);
    try std.testing.expectEqual(3, A[15].im);
    try std.testing.expectEqual(8, A[16].re);
    try std.testing.expectEqual(-1, A[16].im);
    try std.testing.expectEqual(11, A[17].re);
    try std.testing.expectEqual(3, A[17].im);
    try std.testing.expectEqual(239, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(300, A[19].re);
    try std.testing.expectEqual(3, A[19].im);
    try std.testing.expectEqual(5, A[20].re);
    try std.testing.expectEqual(-1, A[20].im);
    try std.testing.expectEqual(9, A[21].re);
    try std.testing.expectEqual(4, A[21].im);
    try std.testing.expectEqual(12, A[22].re);
    try std.testing.expectEqual(-2, A[22].im);
    try std.testing.expectEqual(14, A[23].re);
    try std.testing.expectEqual(1, A[23].im);
    try std.testing.expectEqual(377, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(7, 8),
        Complex(f64).init(5, 6),
        Complex(f64).init(3, 4),
        Complex(f64).init(1, 2),
    });

    blas.her(Complex(f64), .ColumnMajor, .Upper, n, alpha, x2.ptr, -1, A.ptr, n);

    try std.testing.expectEqual(21, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(24, A[1].re);
    try std.testing.expectEqual(3, A[1].im);
    try std.testing.expectEqual(37, A[2].re);
    try std.testing.expectEqual(10, A[2].im);
    try std.testing.expectEqual(50, A[3].re);
    try std.testing.expectEqual(9, A[3].im);
    try std.testing.expectEqual(63, A[4].re);
    try std.testing.expectEqual(17, A[4].im);
    try std.testing.expectEqual(24, A[5].re);
    try std.testing.expectEqual(5, A[5].im);
    try std.testing.expectEqual(106, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(85, A[7].re);
    try std.testing.expectEqual(2, A[7].im);
    try std.testing.expectEqual(114, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(143, A[9].re);
    try std.testing.expectEqual(8, A[9].im);
    try std.testing.expectEqual(37, A[10].re);
    try std.testing.expectEqual(6, A[10].im);
    try std.testing.expectEqual(85, A[11].re);
    try std.testing.expectEqual(6, A[11].im);
    try std.testing.expectEqual(254, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(177, A[13].re);
    try std.testing.expectEqual(1, A[13].im);
    try std.testing.expectEqual(222, A[14].re);
    try std.testing.expectEqual(10, A[14].im);
    try std.testing.expectEqual(50, A[15].re);
    try std.testing.expectEqual(15, A[15].im);
    try std.testing.expectEqual(114, A[16].re);
    try std.testing.expectEqual(7, A[16].im);
    try std.testing.expectEqual(177, A[17].re);
    try std.testing.expectEqual(7, A[17].im);
    try std.testing.expectEqual(465, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(300, A[19].re);
    try std.testing.expectEqual(3, A[19].im);
    try std.testing.expectEqual(63, A[20].re);
    try std.testing.expectEqual(15, A[20].im);
    try std.testing.expectEqual(143, A[21].re);
    try std.testing.expectEqual(16, A[21].im);
    try std.testing.expectEqual(222, A[22].re);
    try std.testing.expectEqual(6, A[22].im);
    try std.testing.expectEqual(300, A[23].re);
    try std.testing.expectEqual(5, A[23].im);
    try std.testing.expectEqual(739, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x3 = try a.alloc(Complex(f64), n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });

    blas.her(Complex(f64), .RowMajor, .Lower, n, alpha, x3.ptr, 1, A.ptr, n);

    try std.testing.expectEqual(31, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(24, A[1].re);
    try std.testing.expectEqual(3, A[1].im);
    try std.testing.expectEqual(37, A[2].re);
    try std.testing.expectEqual(10, A[2].im);
    try std.testing.expectEqual(50, A[3].re);
    try std.testing.expectEqual(9, A[3].im);
    try std.testing.expectEqual(63, A[4].re);
    try std.testing.expectEqual(17, A[4].im);
    try std.testing.expectEqual(46, A[5].re);
    try std.testing.expectEqual(1, A[5].im);
    try std.testing.expectEqual(156, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(85, A[7].re);
    try std.testing.expectEqual(2, A[7].im);
    try std.testing.expectEqual(114, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(143, A[9].re);
    try std.testing.expectEqual(8, A[9].im);
    try std.testing.expectEqual(71, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(163, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(376, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(177, A[13].re);
    try std.testing.expectEqual(1, A[13].im);
    try std.testing.expectEqual(222, A[14].re);
    try std.testing.expectEqual(10, A[14].im);
    try std.testing.expectEqual(96, A[15].re);
    try std.testing.expectEqual(3, A[15].im);
    try std.testing.expectEqual(220, A[16].re);
    try std.testing.expectEqual(-1, A[16].im);
    try std.testing.expectEqual(343, A[17].re);
    try std.testing.expectEqual(3, A[17].im);
    try std.testing.expectEqual(691, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(300, A[19].re);
    try std.testing.expectEqual(3, A[19].im);
    try std.testing.expectEqual(121, A[20].re);
    try std.testing.expectEqual(-1, A[20].im);
    try std.testing.expectEqual(277, A[21].re);
    try std.testing.expectEqual(4, A[21].im);
    try std.testing.expectEqual(432, A[22].re);
    try std.testing.expectEqual(-2, A[22].im);
    try std.testing.expectEqual(586, A[23].re);
    try std.testing.expectEqual(1, A[23].im);
    try std.testing.expectEqual(1101, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x4 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x4);

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

    blas.her(Complex(f64), .ColumnMajor, .Lower, n, alpha, x4.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(41, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(-1, A[1].im);
    try std.testing.expectEqual(71, A[2].re);
    try std.testing.expectEqual(2, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(-3, A[3].im);
    try std.testing.expectEqual(121, A[4].re);
    try std.testing.expectEqual(1, A[4].im);
    try std.testing.expectEqual(46, A[5].re);
    try std.testing.expectEqual(1, A[5].im);
    try std.testing.expectEqual(206, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(163, A[7].re);
    try std.testing.expectEqual(-2, A[7].im);
    try std.testing.expectEqual(220, A[8].re);
    try std.testing.expectEqual(1, A[8].im);
    try std.testing.expectEqual(277, A[9].re);
    try std.testing.expectEqual(-4, A[9].im);
    try std.testing.expectEqual(71, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(163, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(498, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(343, A[13].re);
    try std.testing.expectEqual(-3, A[13].im);
    try std.testing.expectEqual(432, A[14].re);
    try std.testing.expectEqual(2, A[14].im);
    try std.testing.expectEqual(96, A[15].re);
    try std.testing.expectEqual(3, A[15].im);
    try std.testing.expectEqual(220, A[16].re);
    try std.testing.expectEqual(-1, A[16].im);
    try std.testing.expectEqual(343, A[17].re);
    try std.testing.expectEqual(3, A[17].im);
    try std.testing.expectEqual(917, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(586, A[19].re);
    try std.testing.expectEqual(-1, A[19].im);
    try std.testing.expectEqual(121, A[20].re);
    try std.testing.expectEqual(-1, A[20].im);
    try std.testing.expectEqual(277, A[21].re);
    try std.testing.expectEqual(4, A[21].im);
    try std.testing.expectEqual(432, A[22].re);
    try std.testing.expectEqual(-2, A[22].im);
    try std.testing.expectEqual(586, A[23].re);
    try std.testing.expectEqual(1, A[23].im);
    try std.testing.expectEqual(1463, A[24].re);
    try std.testing.expectEqual(0, A[24].im);
}
