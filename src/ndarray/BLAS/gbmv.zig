const std = @import("std");
const core = @import("../../core/core.zig");
const Order = @import("../ndarray.zig").Order;
const Transpose = @import("../ndarray.zig").Transpose;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn gbmv(comptime T: type, order: Order, trans: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (m <= 0 or n <= 0 or kl < 0 or ku < 0 or lda < kl + ku + 1) return;

    var M = m;
    var N = n;
    var KL = kl;
    var KU = ku;
    var TRANS = trans;
    if (order == .RowMajor) {
        M = n;
        N = m;
        KL = ku;
        KU = kl;
        TRANS = if (trans == .NoTrans) .Trans else if (trans == .ConjNoTrans) .ConjTrans else if (trans == .Trans) .NoTrans else .ConjNoTrans;
    }

    if (TRANS == .Trans or TRANS == .ConjTrans) {
        const tmp = M;
        M = N;
        N = tmp;
    }

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.gbmv does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (alpha == 0 and beta == 1) return;

            if (alpha == 0) {
                var iy: isize = if (incy < 0) (-M + 1) * incy else 0;
                const Sty = iy + M * incy;
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

            if (TRANS == .NoTrans or TRANS == .ConjNoTrans) {
                var iy: isize = if (incy < 0) (-M + 1) * incy else 0;
                const Sty = iy + M * incy;
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
                var jx: isize = if (incx < 0) (-N + 1) * incx else 0;
                var ky: isize = 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (M - 1 > j + KL) j + KL else M - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    iy = if (incy < 0) (-M + 1) * incy + ky else ky;
                    while (i <= I1) {
                        y[@intCast(iy)] += t0 * A[@intCast(iaij)];

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    if (j >= KU) ky += incy;

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-M + 1) * incy else 0;
                var kx: isize = 0;
                while (j < M) {
                    var t0: T = 0;
                    const k: isize = KU - j;
                    const I0: isize = if (j - KU > 0) j - KU else 0;
                    const I1: isize = if (N - 1 > j + KL) j + KL else N - 1;

                    var i: isize = I0;
                    var iaij: isize = k + I0 + jaj;
                    var ix: isize = if (incx < 0) (-N + 1) * incx + kx else kx;
                    while (i <= I1) {
                        t0 += A[@intCast(iaij)] * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    if (beta == 0) {
                        y[@intCast(jy)] = 0;
                    } else if (beta != 1) {
                        y[@intCast(jy)] *= beta;
                    }

                    y[@intCast(jy)] += alpha * t0;

                    if (j >= KU) kx += incx;

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            }
        },
        .Complex => {},
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.gbmv only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "gbmv" {
    const a = std.testing.allocator;
    //const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const kl = 1;
    const ku = 2;
    const alpha: f64 = 1;
    const beta: f64 = 1;

    const A = try a.alloc(f64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, m);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]f64{ 0, 0, 3, 7, 11, 0, 2, 6, 10, 14, 1, 5, 9, 13, 0, 4, 8, 12, 0, 0 });
    @memcpy(x1.ptr, &[_]f64{ 1, 2, 3, 4, 5 });
    @memcpy(y1.ptr, &[_]f64{ 1, 2, 3, 4 });

    BLAS.gbmv(f64, .RowMajor, .NoTrans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x1.ptr, 1, beta, y1.ptr, 1);

    try std.testing.expectEqual(28, y1[0]);
    try std.testing.expectEqual(43, y1[1]);
    try std.testing.expectEqual(94, y1[2]);
    try std.testing.expectEqual(83, y1[3]);

    const x2 = try a.alloc(f64, n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * m);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]f64{ 5, 4, 3, 2, 1 });
    @memcpy(y2.ptr, &[_]f64{ 1, 0, 2, 0, 3, 0, 4, 0 });

    BLAS.gbmv(f64, .ColumnMajor, .NoTrans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x2.ptr, -1, beta, y2.ptr, 2);

    try std.testing.expectEqual(34, y2[0]);
    try std.testing.expectEqual(91, y2[2]);
    try std.testing.expectEqual(110, y2[4]);
    try std.testing.expectEqual(79, y2[6]);

    const x3 = try a.alloc(f64, m);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]f64{ 1, 2, 3, 4 });
    @memcpy(y3.ptr, &[_]f64{ 5, 0, 4, 0, 3, 0, 2, 0, 1, 0 });

    BLAS.gbmv(f64, .RowMajor, .Trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x3.ptr, 1, beta, y3.ptr, -2);

    try std.testing.expectEqual(23, y3[8]);
    try std.testing.expectEqual(35, y3[6]);
    try std.testing.expectEqual(92, y3[4]);
    try std.testing.expectEqual(71, y3[2]);
    try std.testing.expectEqual(20, y3[0]);

    const x4 = try a.alloc(f64, 2 * m);
    defer a.free(x4);
    const y4 = try a.alloc(f64, n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]f64{ 4, 0, 3, 0, 2, 0, 1, 0 });
    @memcpy(y4.ptr, &[_]f64{ 5, 4, 3, 2, 1 });

    BLAS.gbmv(f64, .ColumnMajor, .Trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x4.ptr, -2, beta, y4.ptr, -1);

    try std.testing.expectEqual(18, y4[4]);
    try std.testing.expectEqual(24, y4[3]);
    try std.testing.expectEqual(64, y4[2]);
    try std.testing.expectEqual(61, y4[1]);
    try std.testing.expectEqual(77, y4[0]);
}
