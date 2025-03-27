const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;

pub inline fn gerc(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

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

    switch (numericType) {
        .bool => @compileError("blas.gerc does not support bool."),
        .int, .float => @compileError("blas.gerc does not support integer or float types."),
        .cfloat => {
            if (alpha.re == 0 and alpha.im == 0) return;

            if (order == .ColumnMajor) {
                var j: isize = 0;
                var jaj: isize = 0;
                var jy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                while (j < N) {
                    const t0 = T.init(alpha.re * y[@intCast(jy)].re + alpha.im * y[@intCast(jy)].im, alpha.im * y[@intCast(jy)].re - alpha.re * y[@intCast(jy)].im);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var ix: isize = if (incx < 0) (-LENX + 1) * incx else 0;
                    while (i < M) {
                        A[@intCast(iaij)].re += x[@intCast(ix)].re * t0.re - x[@intCast(ix)].im * t0.im;
                        A[@intCast(iaij)].im += x[@intCast(ix)].re * t0.im + x[@intCast(ix)].im * t0.re;

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
                    const t0 = T.init(alpha.re * x[@intCast(jx)].re - alpha.im * x[@intCast(jx)].im, alpha.im * x[@intCast(jx)].re + alpha.re * x[@intCast(jx)].im);

                    var i: isize = 0;
                    var iaij: isize = jaj;
                    var iy: isize = if (incy < 0) (-LENY + 1) * incy else 0;
                    while (i < M) {
                        A[@intCast(iaij)].re += y[@intCast(iy)].re * t0.re + y[@intCast(iy)].im * t0.im;
                        A[@intCast(iaij)].im += y[@intCast(iy)].re * t0.im - y[@intCast(iy)].im * t0.re;

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
        .integer, .rational, .real, .complex, .expression => @compileError("blas.gerc only supports simple types."),
        .unsupported => unreachable,
    }
}

test gerc {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const m = 4;
    const n = 5;
    const alpha: Complex(f64) = Complex(f64).init(1, 1);

    const A = try a.alloc(Complex(f64), m * n);
    defer a.free(A);
    const x1 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x1);
    const y1 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(2, 2),
        Complex(f64).init(3, 3),
        Complex(f64).init(4, 4),
        Complex(f64).init(5, 5),
        Complex(f64).init(6, 6),
        Complex(f64).init(7, 7),
        Complex(f64).init(8, 8),
        Complex(f64).init(9, 9),
        Complex(f64).init(10, 10),
        Complex(f64).init(11, 11),
        Complex(f64).init(12, 12),
        Complex(f64).init(13, 13),
        Complex(f64).init(14, 14),
        Complex(f64).init(15, 15),
        Complex(f64).init(16, 16),
        Complex(f64).init(17, 17),
        Complex(f64).init(18, 18),
        Complex(f64).init(19, 19),
        Complex(f64).init(20, 20),
    });
    @memcpy(x1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
    });
    @memcpy(y1.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    blas.gerc(Complex(f64), .RowMajor, m, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(3, A[0].re);
    try std.testing.expectEqual(5, A[0].im);
    try std.testing.expectEqual(6, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(9, A[2].re);
    try std.testing.expectEqual(15, A[2].im);
    try std.testing.expectEqual(12, A[3].re);
    try std.testing.expectEqual(20, A[3].im);
    try std.testing.expectEqual(15, A[4].re);
    try std.testing.expectEqual(25, A[4].im);
    try std.testing.expectEqual(12, A[5].re);
    try std.testing.expectEqual(14, A[5].im);
    try std.testing.expectEqual(19, A[6].re);
    try std.testing.expectEqual(23, A[6].im);
    try std.testing.expectEqual(26, A[7].re);
    try std.testing.expectEqual(32, A[7].im);
    try std.testing.expectEqual(33, A[8].re);
    try std.testing.expectEqual(41, A[8].im);
    try std.testing.expectEqual(40, A[9].re);
    try std.testing.expectEqual(50, A[9].im);
    try std.testing.expectEqual(21, A[10].re);
    try std.testing.expectEqual(23, A[10].im);
    try std.testing.expectEqual(32, A[11].re);
    try std.testing.expectEqual(36, A[11].im);
    try std.testing.expectEqual(43, A[12].re);
    try std.testing.expectEqual(49, A[12].im);
    try std.testing.expectEqual(54, A[13].re);
    try std.testing.expectEqual(62, A[13].im);
    try std.testing.expectEqual(65, A[14].re);
    try std.testing.expectEqual(75, A[14].im);
    try std.testing.expectEqual(30, A[15].re);
    try std.testing.expectEqual(32, A[15].im);
    try std.testing.expectEqual(45, A[16].re);
    try std.testing.expectEqual(49, A[16].im);
    try std.testing.expectEqual(60, A[17].re);
    try std.testing.expectEqual(66, A[17].im);
    try std.testing.expectEqual(75, A[18].re);
    try std.testing.expectEqual(83, A[18].im);
    try std.testing.expectEqual(90, A[19].re);
    try std.testing.expectEqual(100, A[19].im);

    const x2 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x2);
    const y2 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]Complex(f64){
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
    });
    @memcpy(y2.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
    });

    blas.gerc(Complex(f64), .RowMajor, m, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr, m);

    try std.testing.expectEqual(3, A[0].re);
    try std.testing.expectEqual(5, A[0].im);
    try std.testing.expectEqual(6, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(9, A[2].re);
    try std.testing.expectEqual(15, A[2].im);
    try std.testing.expectEqual(12, A[3].re);
    try std.testing.expectEqual(20, A[3].im);
    try std.testing.expectEqual(15, A[4].re);
    try std.testing.expectEqual(25, A[4].im);
    try std.testing.expectEqual(12, A[5].re);
    try std.testing.expectEqual(14, A[5].im);
    try std.testing.expectEqual(19, A[6].re);
    try std.testing.expectEqual(23, A[6].im);
    try std.testing.expectEqual(26, A[7].re);
    try std.testing.expectEqual(32, A[7].im);
    try std.testing.expectEqual(33, A[8].re);
    try std.testing.expectEqual(41, A[8].im);
    try std.testing.expectEqual(40, A[9].re);
    try std.testing.expectEqual(50, A[9].im);
    try std.testing.expectEqual(21, A[10].re);
    try std.testing.expectEqual(23, A[10].im);
    try std.testing.expectEqual(32, A[11].re);
    try std.testing.expectEqual(36, A[11].im);
    try std.testing.expectEqual(43, A[12].re);
    try std.testing.expectEqual(49, A[12].im);
    try std.testing.expectEqual(54, A[13].re);
    try std.testing.expectEqual(62, A[13].im);
    try std.testing.expectEqual(65, A[14].re);
    try std.testing.expectEqual(75, A[14].im);
    try std.testing.expectEqual(30, A[15].re);
    try std.testing.expectEqual(32, A[15].im);
    try std.testing.expectEqual(45, A[16].re);
    try std.testing.expectEqual(49, A[16].im);
    try std.testing.expectEqual(60, A[17].re);
    try std.testing.expectEqual(66, A[17].im);
    try std.testing.expectEqual(75, A[18].re);
    try std.testing.expectEqual(83, A[18].im);
    try std.testing.expectEqual(90, A[19].re);
    try std.testing.expectEqual(100, A[19].im);

    const x3 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x3);
    const y3 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 6),
        Complex(f64).init(0, 0),
        Complex(f64).init(7, 8),
        Complex(f64).init(0, 0),
    });
    @memcpy(y3.ptr, &[_]Complex(f64){
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
    });

    blas.gerc(Complex(f64), .ColumnMajor, m, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(5, A[0].re);
    try std.testing.expectEqual(9, A[0].im);
    try std.testing.expectEqual(12, A[1].re);
    try std.testing.expectEqual(18, A[1].im);
    try std.testing.expectEqual(19, A[2].re);
    try std.testing.expectEqual(27, A[2].im);
    try std.testing.expectEqual(26, A[3].re);
    try std.testing.expectEqual(36, A[3].im);
    try std.testing.expectEqual(15, A[4].re);
    try std.testing.expectEqual(25, A[4].im);
    try std.testing.expectEqual(16, A[5].re);
    try std.testing.expectEqual(22, A[5].im);
    try std.testing.expectEqual(31, A[6].re);
    try std.testing.expectEqual(39, A[6].im);
    try std.testing.expectEqual(46, A[7].re);
    try std.testing.expectEqual(56, A[7].im);
    try std.testing.expectEqual(61, A[8].re);
    try std.testing.expectEqual(73, A[8].im);
    try std.testing.expectEqual(40, A[9].re);
    try std.testing.expectEqual(50, A[9].im);
    try std.testing.expectEqual(27, A[10].re);
    try std.testing.expectEqual(35, A[10].im);
    try std.testing.expectEqual(50, A[11].re);
    try std.testing.expectEqual(60, A[11].im);
    try std.testing.expectEqual(73, A[12].re);
    try std.testing.expectEqual(85, A[12].im);
    try std.testing.expectEqual(96, A[13].re);
    try std.testing.expectEqual(110, A[13].im);
    try std.testing.expectEqual(65, A[14].re);
    try std.testing.expectEqual(75, A[14].im);
    try std.testing.expectEqual(38, A[15].re);
    try std.testing.expectEqual(48, A[15].im);
    try std.testing.expectEqual(69, A[16].re);
    try std.testing.expectEqual(81, A[16].im);
    try std.testing.expectEqual(100, A[17].re);
    try std.testing.expectEqual(114, A[17].im);
    try std.testing.expectEqual(131, A[18].re);
    try std.testing.expectEqual(147, A[18].im);
    try std.testing.expectEqual(90, A[19].re);
    try std.testing.expectEqual(100, A[19].im);

    const x4 = try a.alloc(Complex(f64), 2 * m);
    defer a.free(x4);
    const y4 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]Complex(f64){
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
        Complex(f64).init(5, 5),
        Complex(f64).init(0, 0),
        Complex(f64).init(4, 4),
        Complex(f64).init(0, 0),
        Complex(f64).init(3, 3),
        Complex(f64).init(0, 0),
        Complex(f64).init(2, 2),
        Complex(f64).init(0, 0),
        Complex(f64).init(1, 1),
        Complex(f64).init(0, 0),
    });

    blas.gerc(Complex(f64), .ColumnMajor, m, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr, m);

    try std.testing.expectEqual(7, A[0].re);
    try std.testing.expectEqual(13, A[0].im);
    try std.testing.expectEqual(18, A[1].re);
    try std.testing.expectEqual(26, A[1].im);
    try std.testing.expectEqual(29, A[2].re);
    try std.testing.expectEqual(39, A[2].im);
    try std.testing.expectEqual(40, A[3].re);
    try std.testing.expectEqual(52, A[3].im);
    try std.testing.expectEqual(19, A[4].re);
    try std.testing.expectEqual(33, A[4].im);
    try std.testing.expectEqual(28, A[5].re);
    try std.testing.expectEqual(38, A[5].im);
    try std.testing.expectEqual(51, A[6].re);
    try std.testing.expectEqual(63, A[6].im);
    try std.testing.expectEqual(74, A[7].re);
    try std.testing.expectEqual(88, A[7].im);
    try std.testing.expectEqual(67, A[8].re);
    try std.testing.expectEqual(85, A[8].im);
    try std.testing.expectEqual(58, A[9].re);
    try std.testing.expectEqual(74, A[9].im);
    try std.testing.expectEqual(57, A[10].re);
    try std.testing.expectEqual(71, A[10].im);
    try std.testing.expectEqual(92, A[11].re);
    try std.testing.expectEqual(108, A[11].im);
    try std.testing.expectEqual(81, A[12].re);
    try std.testing.expectEqual(101, A[12].im);
    try std.testing.expectEqual(120, A[13].re);
    try std.testing.expectEqual(142, A[13].im);
    try std.testing.expectEqual(105, A[14].re);
    try std.testing.expectEqual(123, A[14].im);
    try std.testing.expectEqual(94, A[15].re);
    try std.testing.expectEqual(112, A[15].im);
    try std.testing.expectEqual(79, A[16].re);
    try std.testing.expectEqual(101, A[16].im);
    try std.testing.expectEqual(130, A[17].re);
    try std.testing.expectEqual(154, A[17].im);
    try std.testing.expectEqual(181, A[18].re);
    try std.testing.expectEqual(207, A[18].im);
    try std.testing.expectEqual(160, A[19].re);
    try std.testing.expectEqual(180, A[19].im);
}
