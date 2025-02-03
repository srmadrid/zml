const std = @import("std");
const core = @import("../core/core.zig");
const Order = @import("../ndarray/ndarray.zig").Order;
const Uplo = @import("../ndarray/ndarray.zig").Uplo;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn spr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    const LENX = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.spr does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < j) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    Ap[@intCast(iaij)] += t0 * x[@intCast(jx)];

                    j += 1;
                    jaj += j;
                    jx += incx;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    Ap[@intCast(jaj)] += t0 * x[@intCast(jx)];

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    while (i < I1) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)];

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
        .Complex => @compileError("BLAS.spr does not support complex numbers."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.spr only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "spr" {
    const a = std.testing.allocator;

    const n = 5;
    const alpha: f64 = 2;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);

    @memcpy(A.ptr, &[_]f64{
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    });
    @memcpy(x1.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.spr(f64, .RowMajor, .Upper, n, alpha, x1.ptr, 2, A.ptr);

    try std.testing.expectEqual(3, A[0]);
    try std.testing.expectEqual(6, A[1]);
    try std.testing.expectEqual(9, A[2]);
    try std.testing.expectEqual(12, A[3]);
    try std.testing.expectEqual(15, A[4]);
    try std.testing.expectEqual(14, A[5]);
    try std.testing.expectEqual(19, A[6]);
    try std.testing.expectEqual(24, A[7]);
    try std.testing.expectEqual(29, A[8]);
    try std.testing.expectEqual(28, A[9]);
    try std.testing.expectEqual(35, A[10]);
    try std.testing.expectEqual(42, A[11]);
    try std.testing.expectEqual(45, A[12]);
    try std.testing.expectEqual(54, A[13]);
    try std.testing.expectEqual(65, A[14]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.spr(f64, .RowMajor, .Upper, n, alpha, x2.ptr, -2, A.ptr);

    try std.testing.expectEqual(5, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(22, A[5]);
    try std.testing.expectEqual(31, A[6]);
    try std.testing.expectEqual(40, A[7]);
    try std.testing.expectEqual(49, A[8]);
    try std.testing.expectEqual(46, A[9]);
    try std.testing.expectEqual(59, A[10]);
    try std.testing.expectEqual(72, A[11]);
    try std.testing.expectEqual(77, A[12]);
    try std.testing.expectEqual(94, A[13]);
    try std.testing.expectEqual(115, A[14]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.spr(f64, .ColumnMajor, .Upper, n, alpha, x3.ptr, 2, A.ptr);

    try std.testing.expectEqual(7, A[0]);
    try std.testing.expectEqual(14, A[1]);
    try std.testing.expectEqual(23, A[2]);
    try std.testing.expectEqual(26, A[3]);
    try std.testing.expectEqual(37, A[4]);
    try std.testing.expectEqual(40, A[5]);
    try std.testing.expectEqual(39, A[6]);
    try std.testing.expectEqual(56, A[7]);
    try std.testing.expectEqual(73, A[8]);
    try std.testing.expectEqual(78, A[9]);
    try std.testing.expectEqual(69, A[10]);
    try std.testing.expectEqual(92, A[11]);
    try std.testing.expectEqual(107, A[12]);
    try std.testing.expectEqual(134, A[13]);
    try std.testing.expectEqual(165, A[14]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);

    @memcpy(x4.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.spr(f64, .ColumnMajor, .Upper, n, alpha, x4.ptr, -2, A.ptr);

    try std.testing.expectEqual(9, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(31, A[2]);
    try std.testing.expectEqual(32, A[3]);
    try std.testing.expectEqual(49, A[4]);
    try std.testing.expectEqual(58, A[5]);
    try std.testing.expectEqual(47, A[6]);
    try std.testing.expectEqual(72, A[7]);
    try std.testing.expectEqual(97, A[8]);
    try std.testing.expectEqual(110, A[9]);
    try std.testing.expectEqual(79, A[10]);
    try std.testing.expectEqual(112, A[11]);
    try std.testing.expectEqual(137, A[12]);
    try std.testing.expectEqual(174, A[13]);
    try std.testing.expectEqual(215, A[14]);

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);

    @memcpy(x5.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.spr(f64, .RowMajor, .Lower, n, alpha, x5.ptr, 2, A.ptr);

    try std.testing.expectEqual(11, A[0]);
    try std.testing.expectEqual(22, A[1]);
    try std.testing.expectEqual(39, A[2]);
    try std.testing.expectEqual(38, A[3]);
    try std.testing.expectEqual(61, A[4]);
    try std.testing.expectEqual(76, A[5]);
    try std.testing.expectEqual(55, A[6]);
    try std.testing.expectEqual(88, A[7]);
    try std.testing.expectEqual(121, A[8]);
    try std.testing.expectEqual(142, A[9]);
    try std.testing.expectEqual(89, A[10]);
    try std.testing.expectEqual(132, A[11]);
    try std.testing.expectEqual(167, A[12]);
    try std.testing.expectEqual(214, A[13]);
    try std.testing.expectEqual(265, A[14]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);

    @memcpy(x6.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.spr(f64, .RowMajor, .Lower, n, alpha, x6.ptr, -2, A.ptr);

    try std.testing.expectEqual(13, A[0]);
    try std.testing.expectEqual(26, A[1]);
    try std.testing.expectEqual(47, A[2]);
    try std.testing.expectEqual(44, A[3]);
    try std.testing.expectEqual(73, A[4]);
    try std.testing.expectEqual(94, A[5]);
    try std.testing.expectEqual(63, A[6]);
    try std.testing.expectEqual(104, A[7]);
    try std.testing.expectEqual(145, A[8]);
    try std.testing.expectEqual(174, A[9]);
    try std.testing.expectEqual(99, A[10]);
    try std.testing.expectEqual(152, A[11]);
    try std.testing.expectEqual(197, A[12]);
    try std.testing.expectEqual(254, A[13]);
    try std.testing.expectEqual(315, A[14]);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);

    @memcpy(x7.ptr, &[_]f64{
        1,
        0,
        2,
        0,
        3,
        0,
        4,
        0,
        5,
        0,
    });

    BLAS.spr(f64, .ColumnMajor, .Upper, n, alpha, x7.ptr, 2, A.ptr);

    try std.testing.expectEqual(15, A[0]);
    try std.testing.expectEqual(30, A[1]);
    try std.testing.expectEqual(55, A[2]);
    try std.testing.expectEqual(50, A[3]);
    try std.testing.expectEqual(85, A[4]);
    try std.testing.expectEqual(112, A[5]);
    try std.testing.expectEqual(71, A[6]);
    try std.testing.expectEqual(120, A[7]);
    try std.testing.expectEqual(169, A[8]);
    try std.testing.expectEqual(206, A[9]);
    try std.testing.expectEqual(109, A[10]);
    try std.testing.expectEqual(172, A[11]);
    try std.testing.expectEqual(227, A[12]);
    try std.testing.expectEqual(294, A[13]);
    try std.testing.expectEqual(365, A[14]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);

    @memcpy(x8.ptr, &[_]f64{
        5,
        0,
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    BLAS.spr(f64, .ColumnMajor, .Upper, n, alpha, x8.ptr, -2, A.ptr);

    try std.testing.expectEqual(17, A[0]);
    try std.testing.expectEqual(34, A[1]);
    try std.testing.expectEqual(63, A[2]);
    try std.testing.expectEqual(56, A[3]);
    try std.testing.expectEqual(97, A[4]);
    try std.testing.expectEqual(130, A[5]);
    try std.testing.expectEqual(79, A[6]);
    try std.testing.expectEqual(136, A[7]);
    try std.testing.expectEqual(193, A[8]);
    try std.testing.expectEqual(238, A[9]);
    try std.testing.expectEqual(119, A[10]);
    try std.testing.expectEqual(192, A[11]);
    try std.testing.expectEqual(257, A[12]);
    try std.testing.expectEqual(334, A[13]);
    try std.testing.expectEqual(415, A[14]);
}
