const std = @import("std");
const core = @import("../../core.zig");
const blas = @import("../blas.zig");
const Order = blas.Order;
const Uplo = blas.Uplo;

const Numeric = core.types.Numeric;

pub inline fn her(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Numeric(T), x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    @setRuntimeSafety(false);
    const numericType = core.types.numericType(T);

    if (n <= 0) return;

    const N = n;
    var UPLO = uplo;
    var conj: Numeric(T) = 1;
    if (order == .RowMajor) {
        UPLO = if (uplo == .Upper) .Lower else .Upper;
        conj = -1;
    }

    if (lda < @max(1, N)) return;

    const LENX = N;

    switch (numericType) {
        .bool => @compileError("blas.her does not support bool."),
        .int, .float => @compileError("blas.her does not support int or float."),
        .cfloat => {
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
        .integer, .rational, .real, .complex, .expression => @compileError("blas.her only supports simple types."),
        .unsupported => unreachable,
    }
}

test "her" {
    const a = std.testing.allocator;
    const Complex = std.math.Complex;

    const n = 5;
    const alpha: f64 = 2;

    const A = try a.alloc(Complex(f64), n * n);
    defer a.free(A);
    const x1 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x1);

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
        Complex(f64).init(21, 21),
        Complex(f64).init(22, 22),
        Complex(f64).init(23, 23),
        Complex(f64).init(24, 24),
        Complex(f64).init(25, 25),
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
        Complex(f64).init(9, 10),
        Complex(f64).init(0, 0),
    });

    blas.her(Complex(f64), .RowMajor, .Upper, n, alpha, x1.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(11, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(24, A[1].re);
    try std.testing.expectEqual(6, A[1].im);
    try std.testing.expectEqual(37, A[2].re);
    try std.testing.expectEqual(11, A[2].im);
    try std.testing.expectEqual(50, A[3].re);
    try std.testing.expectEqual(16, A[3].im);
    try std.testing.expectEqual(63, A[4].re);
    try std.testing.expectEqual(21, A[4].im);
    try std.testing.expectEqual(6, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(57, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(86, A[7].re);
    try std.testing.expectEqual(12, A[7].im);
    try std.testing.expectEqual(115, A[8].re);
    try std.testing.expectEqual(17, A[8].im);
    try std.testing.expectEqual(144, A[9].re);
    try std.testing.expectEqual(22, A[9].im);
    try std.testing.expectEqual(11, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(12, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(135, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(180, A[13].re);
    try std.testing.expectEqual(18, A[13].im);
    try std.testing.expectEqual(225, A[14].re);
    try std.testing.expectEqual(23, A[14].im);
    try std.testing.expectEqual(16, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(17, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(18, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(245, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(306, A[19].re);
    try std.testing.expectEqual(24, A[19].im);
    try std.testing.expectEqual(21, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(22, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(23, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(24, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(387, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x2 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x2);

    @memcpy(x2.ptr, &[_]Complex(f64){
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

    blas.her(Complex(f64), .RowMajor, .Upper, n, alpha, x2.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(21, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(71, A[2].re);
    try std.testing.expectEqual(19, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(28, A[3].im);
    try std.testing.expectEqual(121, A[4].re);
    try std.testing.expectEqual(37, A[4].im);
    try std.testing.expectEqual(6, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(107, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(164, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(221, A[8].re);
    try std.testing.expectEqual(25, A[8].im);
    try std.testing.expectEqual(278, A[9].re);
    try std.testing.expectEqual(34, A[9].im);
    try std.testing.expectEqual(11, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(12, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(257, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(346, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(435, A[14].re);
    try std.testing.expectEqual(31, A[14].im);
    try std.testing.expectEqual(16, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(17, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(18, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(471, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(592, A[19].re);
    try std.testing.expectEqual(28, A[19].im);
    try std.testing.expectEqual(21, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(22, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(23, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(24, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(749, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x3 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x3);

    @memcpy(x3.ptr, &[_]Complex(f64){
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

    blas.her(Complex(f64), .ColumnMajor, .Upper, n, alpha, x3.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(31, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(71, A[2].re);
    try std.testing.expectEqual(19, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(28, A[3].im);
    try std.testing.expectEqual(121, A[4].re);
    try std.testing.expectEqual(37, A[4].im);
    try std.testing.expectEqual(28, A[5].re);
    try std.testing.expectEqual(10, A[5].im);
    try std.testing.expectEqual(157, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(164, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(221, A[8].re);
    try std.testing.expectEqual(25, A[8].im);
    try std.testing.expectEqual(278, A[9].re);
    try std.testing.expectEqual(34, A[9].im);
    try std.testing.expectEqual(45, A[10].re);
    try std.testing.expectEqual(19, A[10].im);
    try std.testing.expectEqual(90, A[11].re);
    try std.testing.expectEqual(16, A[11].im);
    try std.testing.expectEqual(379, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(346, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(435, A[14].re);
    try std.testing.expectEqual(31, A[14].im);
    try std.testing.expectEqual(62, A[15].re);
    try std.testing.expectEqual(28, A[15].im);
    try std.testing.expectEqual(123, A[16].re);
    try std.testing.expectEqual(25, A[16].im);
    try std.testing.expectEqual(184, A[17].re);
    try std.testing.expectEqual(22, A[17].im);
    try std.testing.expectEqual(697, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(592, A[19].re);
    try std.testing.expectEqual(28, A[19].im);
    try std.testing.expectEqual(79, A[20].re);
    try std.testing.expectEqual(37, A[20].im);
    try std.testing.expectEqual(156, A[21].re);
    try std.testing.expectEqual(34, A[21].im);
    try std.testing.expectEqual(233, A[22].re);
    try std.testing.expectEqual(31, A[22].im);
    try std.testing.expectEqual(310, A[23].re);
    try std.testing.expectEqual(28, A[23].im);
    try std.testing.expectEqual(1111, A[24].re);
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

    blas.her(Complex(f64), .ColumnMajor, .Upper, n, alpha, x4.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(41, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(71, A[2].re);
    try std.testing.expectEqual(19, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(28, A[3].im);
    try std.testing.expectEqual(121, A[4].re);
    try std.testing.expectEqual(37, A[4].im);
    try std.testing.expectEqual(50, A[5].re);
    try std.testing.expectEqual(14, A[5].im);
    try std.testing.expectEqual(207, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(164, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(221, A[8].re);
    try std.testing.expectEqual(25, A[8].im);
    try std.testing.expectEqual(278, A[9].re);
    try std.testing.expectEqual(34, A[9].im);
    try std.testing.expectEqual(79, A[10].re);
    try std.testing.expectEqual(27, A[10].im);
    try std.testing.expectEqual(168, A[11].re);
    try std.testing.expectEqual(20, A[11].im);
    try std.testing.expectEqual(501, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(346, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(435, A[14].re);
    try std.testing.expectEqual(31, A[14].im);
    try std.testing.expectEqual(108, A[15].re);
    try std.testing.expectEqual(40, A[15].im);
    try std.testing.expectEqual(229, A[16].re);
    try std.testing.expectEqual(33, A[16].im);
    try std.testing.expectEqual(350, A[17].re);
    try std.testing.expectEqual(26, A[17].im);
    try std.testing.expectEqual(923, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(592, A[19].re);
    try std.testing.expectEqual(28, A[19].im);
    try std.testing.expectEqual(137, A[20].re);
    try std.testing.expectEqual(53, A[20].im);
    try std.testing.expectEqual(290, A[21].re);
    try std.testing.expectEqual(46, A[21].im);
    try std.testing.expectEqual(443, A[22].re);
    try std.testing.expectEqual(39, A[22].im);
    try std.testing.expectEqual(596, A[23].re);
    try std.testing.expectEqual(32, A[23].im);
    try std.testing.expectEqual(1473, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x5 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x5);

    @memcpy(x5.ptr, &[_]Complex(f64){
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

    blas.her(Complex(f64), .RowMajor, .Lower, n, alpha, x5.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(51, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(71, A[2].re);
    try std.testing.expectEqual(19, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(28, A[3].im);
    try std.testing.expectEqual(121, A[4].re);
    try std.testing.expectEqual(37, A[4].im);
    try std.testing.expectEqual(72, A[5].re);
    try std.testing.expectEqual(10, A[5].im);
    try std.testing.expectEqual(257, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(164, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(221, A[8].re);
    try std.testing.expectEqual(25, A[8].im);
    try std.testing.expectEqual(278, A[9].re);
    try std.testing.expectEqual(34, A[9].im);
    try std.testing.expectEqual(113, A[10].re);
    try std.testing.expectEqual(19, A[10].im);
    try std.testing.expectEqual(246, A[11].re);
    try std.testing.expectEqual(16, A[11].im);
    try std.testing.expectEqual(623, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(346, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(435, A[14].re);
    try std.testing.expectEqual(31, A[14].im);
    try std.testing.expectEqual(154, A[15].re);
    try std.testing.expectEqual(28, A[15].im);
    try std.testing.expectEqual(335, A[16].re);
    try std.testing.expectEqual(25, A[16].im);
    try std.testing.expectEqual(516, A[17].re);
    try std.testing.expectEqual(22, A[17].im);
    try std.testing.expectEqual(1149, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(592, A[19].re);
    try std.testing.expectEqual(28, A[19].im);
    try std.testing.expectEqual(195, A[20].re);
    try std.testing.expectEqual(37, A[20].im);
    try std.testing.expectEqual(424, A[21].re);
    try std.testing.expectEqual(34, A[21].im);
    try std.testing.expectEqual(653, A[22].re);
    try std.testing.expectEqual(31, A[22].im);
    try std.testing.expectEqual(882, A[23].re);
    try std.testing.expectEqual(28, A[23].im);
    try std.testing.expectEqual(1835, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x6 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x6);

    @memcpy(x6.ptr, &[_]Complex(f64){
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

    blas.her(Complex(f64), .RowMajor, .Lower, n, alpha, x6.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(61, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(46, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(71, A[2].re);
    try std.testing.expectEqual(19, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(28, A[3].im);
    try std.testing.expectEqual(121, A[4].re);
    try std.testing.expectEqual(37, A[4].im);
    try std.testing.expectEqual(94, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(307, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(164, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(221, A[8].re);
    try std.testing.expectEqual(25, A[8].im);
    try std.testing.expectEqual(278, A[9].re);
    try std.testing.expectEqual(34, A[9].im);
    try std.testing.expectEqual(147, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(324, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(745, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(346, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(435, A[14].re);
    try std.testing.expectEqual(31, A[14].im);
    try std.testing.expectEqual(200, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(441, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(682, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(1375, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(592, A[19].re);
    try std.testing.expectEqual(28, A[19].im);
    try std.testing.expectEqual(253, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(558, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(863, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(1168, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(2197, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x7 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x7);

    @memcpy(x7.ptr, &[_]Complex(f64){
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

    blas.her(Complex(f64), .ColumnMajor, .Lower, n, alpha, x7.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(71, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(68, A[1].re);
    try std.testing.expectEqual(6, A[1].im);
    try std.testing.expectEqual(105, A[2].re);
    try std.testing.expectEqual(11, A[2].im);
    try std.testing.expectEqual(142, A[3].re);
    try std.testing.expectEqual(16, A[3].im);
    try std.testing.expectEqual(179, A[4].re);
    try std.testing.expectEqual(21, A[4].im);
    try std.testing.expectEqual(94, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(357, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(242, A[7].re);
    try std.testing.expectEqual(12, A[7].im);
    try std.testing.expectEqual(327, A[8].re);
    try std.testing.expectEqual(17, A[8].im);
    try std.testing.expectEqual(412, A[9].re);
    try std.testing.expectEqual(22, A[9].im);
    try std.testing.expectEqual(147, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(324, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(867, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(512, A[13].re);
    try std.testing.expectEqual(18, A[13].im);
    try std.testing.expectEqual(645, A[14].re);
    try std.testing.expectEqual(23, A[14].im);
    try std.testing.expectEqual(200, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(441, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(682, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(1601, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(878, A[19].re);
    try std.testing.expectEqual(24, A[19].im);
    try std.testing.expectEqual(253, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(558, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(863, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(1168, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(2559, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

    const x8 = try a.alloc(Complex(f64), 2 * n);
    defer a.free(x8);

    @memcpy(x8.ptr, &[_]Complex(f64){
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

    blas.her(Complex(f64), .ColumnMajor, .Lower, n, alpha, x8.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(81, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(90, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(139, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(188, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(237, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(94, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(407, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(320, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(433, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(546, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(147, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(324, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(989, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(678, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(855, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(200, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(441, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(682, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(1827, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(1164, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(253, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(558, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(863, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(1168, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(2921, A[24].re);
    try std.testing.expectEqual(0, A[24].im);
}
