const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Uplo = @import("../ndarray.zig").Uplo;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn symv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;
    const LENY = N;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.symv does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (alpha == 0 and beta == 1) return;

            if (alpha == 0) {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)] = 0;

                        iy += incy;
                    }
                } else if (beta != 1) {
                    while (iy != Sty) {
                        y[@intCast(iy)] *= beta;

                        iy += incy;
                    }
                }

                return;
            }

            if (UPLO == .Upper) {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)] = 0;

                        iy += incy;
                    }
                } else if (beta != 1) {
                    while (iy != Sty) {
                        y[@intCast(iy)] *= beta;

                        iy += incy;
                    }
                }

                var j: isize = 0;
                var jaj: isize = 0;
                var jx = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    var t1: T = 0;

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix = if (incx < 0) (-LENX + 1) * incx else 0;
                    iy = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < j) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        t1 += A[@intCast(iaij)] * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)] += t0 * A[@intCast(iaij)] + alpha * t1;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                    jy += incy;
                }
            } else {
                var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                const Sty = iy + LENY * incy;
                if (beta == 0) {
                    while (iy != Sty) {
                        y[@intCast(iy)] = 0;

                        iy += incy;
                    }
                } else if (beta != 1) {
                    while (iy != Sty) {
                        y[@intCast(iy)] *= beta;

                        iy += incy;
                    }
                }

                const ldap1 = lda + 1;

                var j: isize = 0;
                var jaj: isize = 0;
                var jx = if (incx < 0) (-LENX + 1) * incx else 0;
                var jy = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    var t1: T = 0;
                    const I0: isize = j + 1;
                    const I1: isize = N;

                    var i: isize = I0;
                    var iaij: isize = jaj + 1;
                    var ix = if (incx < 0) (-LENX + 1) * incx + I0 * incx else I0 * incx;
                    iy = if (incy < 0) (-LENY + 1) * incy + I0 * incy else I0 * incy;
                    while (i < I1) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        t1 += A[@intCast(iaij)] * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                        iy += incy;
                    }

                    y[@intCast(jy)] += t0 * A[@intCast(jaj)] + alpha * t1;

                    j += 1;
                    jaj += ldap1;
                    jx += incx;
                    jy += incy;
                }
            }
        },
        .Complex => @compileError("BLAS.symv does not support complex numbers."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.symv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "symv" {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = 2;
    const beta = 3;

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

    BLAS.symv(f64, .RowMajor, .Upper, n, alpha, A.ptr, n, x1.ptr, 1, beta, y1.ptr, 1);

    try std.testing.expectEqual(113, y1[0]);
    try std.testing.expectEqual(258, y1[1]);
    try std.testing.expectEqual(387, y1[2]);
    try std.testing.expectEqual(492, y1[3]);
    try std.testing.expectEqual(565, y1[4]);

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

    BLAS.symv(f64, .ColumnMajor, .Upper, n, alpha, A.ptr, n, x2.ptr, -1, beta, y2.ptr, 2);

    try std.testing.expectEqual(433, y2[0]);
    try std.testing.expectEqual(474, y2[2]);
    try std.testing.expectEqual(531, y2[4]);
    try std.testing.expectEqual(612, y2[6]);
    try std.testing.expectEqual(725, y2[8]);

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

    BLAS.symv(f64, .RowMajor, .Lower, n, alpha, A.ptr, n, x3.ptr, 1, beta, y3.ptr, -2);

    try std.testing.expectEqual(433, y3[8]);
    try std.testing.expectEqual(474, y3[6]);
    try std.testing.expectEqual(531, y3[4]);
    try std.testing.expectEqual(612, y3[2]);
    try std.testing.expectEqual(725, y3[0]);

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

    BLAS.symv(f64, .ColumnMajor, .Lower, n, alpha, A.ptr, n, x4.ptr, -2, beta, y4.ptr, -1);

    try std.testing.expectEqual(113, y4[4]);
    try std.testing.expectEqual(258, y4[3]);
    try std.testing.expectEqual(387, y4[2]);
    try std.testing.expectEqual(492, y4[1]);
    try std.testing.expectEqual(565, y4[0]);
}