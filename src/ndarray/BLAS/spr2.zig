const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Uplo = @import("../ndarray.zig").Uplo;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn spr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, Ap: [*]T) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    const LENX = N;
    const LENY = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.spr2 does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (alpha == 0) return;

            if (UPLO == .Upper) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const t1 = alpha * y[@intCast(jy)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix = if (incx < 0) (-LENX + 1) * incx else 0;
                    var iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)] + t1 * y[@intCast(iy)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    Ap[@intCast(iaij)] += t0 * x[@intCast(jx)] + t1 * y[@intCast(jy)];

                    j += 1;
                    jaj += j;
                    jx += incx;
                    jy += incy;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const t1 = alpha * y[@intCast(jy)];
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    Ap[@intCast(jaj)] += t0 * x[@intCast(jx)] + t1 * y[@intCast(jy)];

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    var iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        Ap[@intCast(iaij)] += t0 * x[@intCast(ix)] + t1 * y[@intCast(iy)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    j += 1;
                    jaj += N - j + 1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .Complex => @compileError("BLAS.spr2 does not support complex numbers."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.spr2 only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "spr2" {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = 2;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, n);
    defer a.free(y1);

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
        2,
        3,
        4,
        5,
    });
    @memcpy(y1.ptr, &[_]f64{
        1,
        2,
        3,
        4,
        5,
    });

    BLAS.spr2(f64, .RowMajor, .Upper, n, alpha, x1.ptr, 1, y1.ptr, 1, A.ptr);

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

    const x2 = try a.alloc(f64, n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]f64{
        5,
        4,
        3,
        2,
        1,
    });
    @memcpy(y2.ptr, &[_]f64{
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

    BLAS.spr2(f64, .ColumnMajor, .Upper, n, alpha, x2.ptr, -1, y2.ptr, 2, A.ptr);

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

    const x3 = try a.alloc(f64, n);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]f64{
        1,
        2,
        3,
        4,
        5,
    });
    @memcpy(y3.ptr, &[_]f64{
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

    BLAS.spr2(f64, .RowMajor, .Lower, n, alpha, x3.ptr, 1, y3.ptr, -2, A.ptr);

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

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(f64, n);
    defer a.free(y4);

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
    @memcpy(y4.ptr, &[_]f64{
        5,
        4,
        3,
        2,
        1,
    });

    BLAS.spr2(f64, .ColumnMajor, .Lower, n, alpha, x4.ptr, -2, y4.ptr, -1, A.ptr);

    try std.testing.expectEqual(17, A[0]);
    try std.testing.expectEqual(34, A[1]);
    try std.testing.expectEqual(59, A[2]);
    try std.testing.expectEqual(60, A[3]);
    try std.testing.expectEqual(93, A[4]);
    try std.testing.expectEqual(110, A[5]);
    try std.testing.expectEqual(87, A[6]);
    try std.testing.expectEqual(136, A[7]);
    try std.testing.expectEqual(185, A[8]);
    try std.testing.expectEqual(210, A[9]);
    try std.testing.expectEqual(147, A[10]);
    try std.testing.expectEqual(212, A[11]);
    try std.testing.expectEqual(261, A[12]);
    try std.testing.expectEqual(334, A[13]);
    try std.testing.expectEqual(415, A[14]);
}
