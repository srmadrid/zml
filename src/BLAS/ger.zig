const std = @import("std");
const core = @import("../core/core.zig");
const Order = @import("../ndarray/ndarray.zig").Order;
const BLAS = @import("BLAS.zig");

const scalar = core.supported.scalar;

pub inline fn ger(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const supported = core.supported.whatSupportedNumericType(T);

    if (m <= 0 or n <= 0) return;

    var M = m;
    var N = n;
    if (order == .RowMajor) {
        M = n;
        N = m;
    }

    if (lda < @max(1, M)) return;

    const LENX = m;
    const LENY = n;

    switch (supported) {
        .BuiltinBool => @compileError("BLAS.ger does not support bool."),
        .BuiltinInt, .BuiltinFloat => {
            if (alpha == 0) return;

            if (order == .ColumnMajor) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = alpha * y[@intCast(jy)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < M) {
                        A[@intCast(iaij)] += t0 * x[@intCast(ix)];

                        i += 1;
                        iaij += 1;
                        ix += incx;
                    }

                    j += 1;
                    jaj += lda;
                    jy += incy;
                }
            } else {
                var j: isize = 0;
                var jaj: isize = 0;
                var jx: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                while (j < N) {
                    const t0 = alpha * x[@intCast(jx)];

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < M) {
                        A[@intCast(iaij)] += t0 * y[@intCast(iy)];

                        i += 1;
                        iaij += 1;
                        iy += incy;
                    }

                    j += 1;
                    jaj += lda;
                    jx += incx;
                }
            }
        },
        .Complex => @compileError("BLAS.ger does not support complex types."),
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("BLAS.ger only supports simple types."),
        .Unsupported => unreachable,
    }
}

test "ger" {
    const a = std.testing.allocator;

    const m = 4;
    const n = 5;
    const alpha: f64 = 1;

    const A = try a.alloc(f64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, m);
    defer a.free(x1);
    const y1 = try a.alloc(f64, n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]f64{ 1, 0, 3, 7, 11, 4, 2, 6, 10, 14, 1, 5, 9, 13, 10, 4, 8, 12, 1, 3 });
    @memcpy(x1.ptr, &[_]f64{ 1, 2, 3, 4 });
    @memcpy(y1.ptr, &[_]f64{ 1, 2, 3, 4, 5 });

    BLAS.ger(f64, .RowMajor, m, n, alpha, x1.ptr, 1, y1.ptr, 1, A.ptr, n);

    try std.testing.expectEqual(2, A[0]);
    try std.testing.expectEqual(2, A[1]);
    try std.testing.expectEqual(6, A[2]);
    try std.testing.expectEqual(11, A[3]);
    try std.testing.expectEqual(16, A[4]);
    try std.testing.expectEqual(6, A[5]);
    try std.testing.expectEqual(6, A[6]);
    try std.testing.expectEqual(12, A[7]);
    try std.testing.expectEqual(18, A[8]);
    try std.testing.expectEqual(24, A[9]);
    try std.testing.expectEqual(4, A[10]);
    try std.testing.expectEqual(11, A[11]);
    try std.testing.expectEqual(18, A[12]);
    try std.testing.expectEqual(25, A[13]);
    try std.testing.expectEqual(25, A[14]);
    try std.testing.expectEqual(8, A[15]);
    try std.testing.expectEqual(16, A[16]);
    try std.testing.expectEqual(24, A[17]);
    try std.testing.expectEqual(17, A[18]);
    try std.testing.expectEqual(23, A[19]);

    const x2 = try a.alloc(f64, m);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]f64{ 4, 3, 2, 1 });
    @memcpy(y2.ptr, &[_]f64{ 1, 0, 2, 0, 3, 0, 4, 0, 5, 0 });

    BLAS.ger(f64, .ColumnMajor, m, n, alpha, x2.ptr, -1, y2.ptr, 2, A.ptr, m);

    try std.testing.expectEqual(3, A[0]);
    try std.testing.expectEqual(4, A[1]);
    try std.testing.expectEqual(9, A[2]);
    try std.testing.expectEqual(15, A[3]);
    try std.testing.expectEqual(18, A[4]);
    try std.testing.expectEqual(10, A[5]);
    try std.testing.expectEqual(12, A[6]);
    try std.testing.expectEqual(20, A[7]);
    try std.testing.expectEqual(21, A[8]);
    try std.testing.expectEqual(30, A[9]);
    try std.testing.expectEqual(13, A[10]);
    try std.testing.expectEqual(23, A[11]);
    try std.testing.expectEqual(22, A[12]);
    try std.testing.expectEqual(33, A[13]);
    try std.testing.expectEqual(37, A[14]);
    try std.testing.expectEqual(24, A[15]);
    try std.testing.expectEqual(21, A[16]);
    try std.testing.expectEqual(34, A[17]);
    try std.testing.expectEqual(32, A[18]);
    try std.testing.expectEqual(43, A[19]);

    const x3 = try a.alloc(f64, m);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]f64{ 1, 2, 3, 4 });
    @memcpy(y3.ptr, &[_]f64{ 5, 0, 4, 0, 3, 0, 2, 0, 1, 0 });

    BLAS.ger(f64, .RowMajor, m, n, alpha, x3.ptr, 1, y3.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(4, A[0]);
    try std.testing.expectEqual(6, A[1]);
    try std.testing.expectEqual(12, A[2]);
    try std.testing.expectEqual(19, A[3]);
    try std.testing.expectEqual(23, A[4]);
    try std.testing.expectEqual(12, A[5]);
    try std.testing.expectEqual(16, A[6]);
    try std.testing.expectEqual(26, A[7]);
    try std.testing.expectEqual(29, A[8]);
    try std.testing.expectEqual(40, A[9]);
    try std.testing.expectEqual(16, A[10]);
    try std.testing.expectEqual(29, A[11]);
    try std.testing.expectEqual(31, A[12]);
    try std.testing.expectEqual(45, A[13]);
    try std.testing.expectEqual(52, A[14]);
    try std.testing.expectEqual(28, A[15]);
    try std.testing.expectEqual(29, A[16]);
    try std.testing.expectEqual(46, A[17]);
    try std.testing.expectEqual(48, A[18]);
    try std.testing.expectEqual(63, A[19]);

    const x4 = try a.alloc(f64, 2 * m);
    defer a.free(x4);
    const y4 = try a.alloc(f64, n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]f64{ 4, 0, 3, 0, 2, 0, 1, 0 });
    @memcpy(y4.ptr, &[_]f64{ 5, 4, 3, 2, 1 });

    BLAS.ger(f64, .ColumnMajor, m, n, alpha, x4.ptr, -2, y4.ptr, -1, A.ptr, m);

    try std.testing.expectEqual(5, A[0]);
    try std.testing.expectEqual(8, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(23, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(16, A[5]);
    try std.testing.expectEqual(22, A[6]);
    try std.testing.expectEqual(34, A[7]);
    try std.testing.expectEqual(32, A[8]);
    try std.testing.expectEqual(46, A[9]);
    try std.testing.expectEqual(25, A[10]);
    try std.testing.expectEqual(41, A[11]);
    try std.testing.expectEqual(35, A[12]);
    try std.testing.expectEqual(53, A[13]);
    try std.testing.expectEqual(64, A[14]);
    try std.testing.expectEqual(44, A[15]);
    try std.testing.expectEqual(34, A[16]);
    try std.testing.expectEqual(56, A[17]);
    try std.testing.expectEqual(63, A[18]);
    try std.testing.expectEqual(83, A[19]);
}
