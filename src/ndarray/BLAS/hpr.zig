const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Uplo = @import("../ndarray.zig").Uplo;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn hpr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
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

    const LENX = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.hpr does not support bool."),
        .BuiltinInt, .BuiltinFloat => @compileError("BLAS.hpr does not support int or float."),
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
                        Ap[@intCast(iaij)].re += t0.re * x[@intCast(ix)].re - t0.im * x[@intCast(ix)].im * conj;
                        Ap[@intCast(iaij)].im += t0.im * x[@intCast(ix)].re + t0.re * x[@intCast(ix)].im * conj;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    Ap[@intCast(iaij)].re += t0.re * x[@intCast(jx)].re - t0.im * x[@intCast(jx)].im * conj;
                    Ap[@intCast(iaij)].im = 0;

                    j += 1;
                    jaj += j;
                    jx += incx;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = T.init(alpha * x[@intCast(jx)].re, alpha * x[@intCast(jx)].im * (-conj));
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    Ap[@intCast(jaj)].re += t0.re * x[@intCast(jx)].re - t0.im * x[@intCast(jx)].im * conj;
                    Ap[@intCast(jaj)].im = 0;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    while (i < I1) {
                        Ap[@intCast(iaij)].re += t0.re * x[@intCast(ix)].re - t0.im * x[@intCast(ix)].im * conj;
                        Ap[@intCast(iaij)].im += t0.im * x[@intCast(ix)].re + t0.re * x[@intCast(ix)].im * conj;

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += N - j + 1;
                    jx += incx;
                }
            }
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.hpr only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "hpr" {
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

    BLAS.hpr(Complex(f64), .RowMajor, .Upper, n, alpha, x1.ptr, 1, A.ptr);

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
    try std.testing.expectEqual(52, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(84, A[6].re);
    try std.testing.expectEqual(4, A[6].im);
    try std.testing.expectEqual(113, A[7].re);
    try std.testing.expectEqual(6, A[7].im);
    try std.testing.expectEqual(142, A[8].re);
    try std.testing.expectEqual(13, A[8].im);
    try std.testing.expectEqual(131, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(169, A[10].re);
    try std.testing.expectEqual(2, A[10].im);
    try std.testing.expectEqual(217, A[11].re);
    try std.testing.expectEqual(10, A[11].im);
    try std.testing.expectEqual(236, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(297, A[13].re);
    try std.testing.expectEqual(1, A[13].im);
    try std.testing.expectEqual(374, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

    const x2 = try a.alloc(Complex(f64), n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]Complex(f64){
        Complex(f64).init(9, 10),
        Complex(f64).init(7, 8),
        Complex(f64).init(5, 6),
        Complex(f64).init(3, 4),
        Complex(f64).init(1, 2),
    });

    BLAS.hpr(Complex(f64), .ColumnMajor, .Upper, n, alpha, x2.ptr, -1, A.ptr);

    try std.testing.expectEqual(21, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(7, A[1].im);
    try std.testing.expectEqual(87, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(84, A[3].re);
    try std.testing.expectEqual(17, A[3].im);
    try std.testing.expectEqual(141, A[4].re);
    try std.testing.expectEqual(21, A[4].im);
    try std.testing.expectEqual(174, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(130, A[6].re);
    try std.testing.expectEqual(16, A[6].im);
    try std.testing.expectEqual(219, A[7].re);
    try std.testing.expectEqual(14, A[7].im);
    try std.testing.expectEqual(308, A[8].re);
    try std.testing.expectEqual(17, A[8].im);
    try std.testing.expectEqual(357, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(227, A[10].re);
    try std.testing.expectEqual(18, A[10].im);
    try std.testing.expectEqual(351, A[11].re);
    try std.testing.expectEqual(22, A[11].im);
    try std.testing.expectEqual(446, A[12].re);
    try std.testing.expectEqual(8, A[12].im);
    try std.testing.expectEqual(583, A[13].re);
    try std.testing.expectEqual(5, A[13].im);
    try std.testing.expectEqual(736, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

    const x3 = try a.alloc(Complex(f64), n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(3, 4),
        Complex(f64).init(5, 6),
        Complex(f64).init(7, 8),
        Complex(f64).init(9, 10),
    });

    BLAS.hpr(Complex(f64), .RowMajor, .Lower, n, alpha, x3.ptr, 1, A.ptr);

    try std.testing.expectEqual(31, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(68, A[1].re);
    try std.testing.expectEqual(3, A[1].im);
    try std.testing.expectEqual(137, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(118, A[3].re);
    try std.testing.expectEqual(9, A[3].im);
    try std.testing.expectEqual(219, A[4].re);
    try std.testing.expectEqual(17, A[4].im);
    try std.testing.expectEqual(296, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(176, A[6].re);
    try std.testing.expectEqual(4, A[6].im);
    try std.testing.expectEqual(325, A[7].re);
    try std.testing.expectEqual(6, A[7].im);
    try std.testing.expectEqual(474, A[8].re);
    try std.testing.expectEqual(13, A[8].im);
    try std.testing.expectEqual(583, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(285, A[10].re);
    try std.testing.expectEqual(2, A[10].im);
    try std.testing.expectEqual(485, A[11].re);
    try std.testing.expectEqual(10, A[11].im);
    try std.testing.expectEqual(656, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(869, A[13].re);
    try std.testing.expectEqual(1, A[13].im);
    try std.testing.expectEqual(1098, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    BLAS.hpr(Complex(f64), .ColumnMajor, .Lower, n, alpha, x4.ptr, -2, A.ptr);

    try std.testing.expectEqual(41, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(90, A[1].re);
    try std.testing.expectEqual(-1, A[1].im);
    try std.testing.expectEqual(171, A[2].re);
    try std.testing.expectEqual(-8, A[2].im);
    try std.testing.expectEqual(164, A[3].re);
    try std.testing.expectEqual(-3, A[3].im);
    try std.testing.expectEqual(277, A[4].re);
    try std.testing.expectEqual(1, A[4].im);
    try std.testing.expectEqual(346, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(254, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(431, A[7].re);
    try std.testing.expectEqual(-2, A[7].im);
    try std.testing.expectEqual(608, A[8].re);
    try std.testing.expectEqual(1, A[8].im);
    try std.testing.expectEqual(705, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(451, A[10].re);
    try std.testing.expectEqual(-2, A[10].im);
    try std.testing.expectEqual(695, A[11].re);
    try std.testing.expectEqual(2, A[11].im);
    try std.testing.expectEqual(882, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(1155, A[13].re);
    try std.testing.expectEqual(-3, A[13].im);
    try std.testing.expectEqual(1460, A[14].re);
    try std.testing.expectEqual(0, A[14].im);
}
