const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const hpmv = zml.linalg.blas.hpmv;

test hpmv {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = cf64.init(1, 1);
    const beta = cf64.init(3, 3);

    const A = try a.alloc(cf64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(cf64, 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(cf64, 2 * n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(2, 2),
        cf64.init(3, 3),
        cf64.init(4, 4),
        cf64.init(5, 5),
        cf64.init(6, 6),
        cf64.init(7, 7),
        cf64.init(8, 8),
        cf64.init(9, 9),
        cf64.init(10, 10),
        cf64.init(11, 11),
        cf64.init(12, 12),
        cf64.init(13, 13),
        cf64.init(14, 14),
        cf64.init(15, 15),
        cf64.init(16, 16),
        cf64.init(17, 17),
        cf64.init(18, 18),
        cf64.init(19, 19),
        cf64.init(20, 20),
        cf64.init(21, 21),
        cf64.init(22, 22),
        cf64.init(23, 23),
        cf64.init(24, 24),
        cf64.init(25, 25),
    });
    @memcpy(x1.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(9, 10),
        cf64.init(0, 0),
    });
    @memcpy(y1.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(5, 5),
        cf64.init(0, 0),
    });

    hpmv(cf64, .RowMajor, .Upper, n, alpha, A.ptr, x1.ptr, 2, beta, y1.ptr, 2);

    try std.testing.expectEqual(-217, y1[0].re);
    try std.testing.expectEqual(197, y1[0].im);
    try std.testing.expectEqual(-394, y1[2].re);
    try std.testing.expectEqual(406, y1[2].im);
    try std.testing.expectEqual(-378, y1[4].re);
    try std.testing.expectEqual(566, y1[4].im);
    try std.testing.expectEqual(-127, y1[6].re);
    try std.testing.expectEqual(683, y1[6].im);
    try std.testing.expectEqual(365, y1[8].re);
    try std.testing.expectEqual(775, y1[8].im);

    const x2 = try a.alloc(cf64, 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(cf64, 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]cf64{
        cf64.init(9, 10),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y2.ptr, &[_]cf64{
        cf64.init(5, 5),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    hpmv(cf64, .RowMajor, .Upper, n, alpha, A.ptr, x2.ptr, -2, beta, y2.ptr, -2);

    try std.testing.expectEqual(-217, y2[8].re);
    try std.testing.expectEqual(197, y2[8].im);
    try std.testing.expectEqual(-394, y2[6].re);
    try std.testing.expectEqual(406, y2[6].im);
    try std.testing.expectEqual(-378, y2[4].re);
    try std.testing.expectEqual(566, y2[4].im);
    try std.testing.expectEqual(-127, y2[2].re);
    try std.testing.expectEqual(683, y2[2].im);
    try std.testing.expectEqual(365, y2[0].re);
    try std.testing.expectEqual(775, y2[0].im);

    const x3 = try a.alloc(cf64, 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(cf64, 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(9, 10),
        cf64.init(0, 0),
    });
    @memcpy(y3.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(5, 5),
        cf64.init(0, 0),
    });

    hpmv(cf64, .ColumnMajor, .Upper, n, alpha, A.ptr, x3.ptr, 2, beta, y3.ptr, 2);

    try std.testing.expectEqual(-397, y3[0].re);
    try std.testing.expectEqual(357, y3[0].im);
    try std.testing.expectEqual(-427, y3[2].re);
    try std.testing.expectEqual(419, y3[2].im);
    try std.testing.expectEqual(-372, y3[4].re);
    try std.testing.expectEqual(500, y3[4].im);
    try std.testing.expectEqual(-138, y3[6].re);
    try std.testing.expectEqual(626, y3[6].im);
    try std.testing.expectEqual(405, y3[8].re);
    try std.testing.expectEqual(835, y3[8].im);

    const x4 = try a.alloc(cf64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(cf64, 2 * n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]cf64{
        cf64.init(9, 10),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y4.ptr, &[_]cf64{
        cf64.init(5, 5),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    hpmv(cf64, .ColumnMajor, .Upper, n, alpha, A.ptr, x4.ptr, -2, beta, y4.ptr, -2);

    try std.testing.expectEqual(-397, y4[8].re);
    try std.testing.expectEqual(357, y4[8].im);
    try std.testing.expectEqual(-427, y4[6].re);
    try std.testing.expectEqual(419, y4[6].im);
    try std.testing.expectEqual(-372, y4[4].re);
    try std.testing.expectEqual(500, y4[4].im);
    try std.testing.expectEqual(-138, y4[2].re);
    try std.testing.expectEqual(626, y4[2].im);
    try std.testing.expectEqual(405, y4[0].re);
    try std.testing.expectEqual(835, y4[0].im);

    const x5 = try a.alloc(cf64, 2 * n);
    defer a.free(x5);
    const y5 = try a.alloc(cf64, 2 * n);
    defer a.free(y5);

    @memcpy(x5.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(9, 10),
        cf64.init(0, 0),
    });
    @memcpy(y5.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(5, 5),
        cf64.init(0, 0),
    });

    hpmv(cf64, .RowMajor, .Lower, n, alpha, A.ptr, x5.ptr, 2, beta, y5.ptr, 2);

    try std.testing.expectEqual(347, y5[0].re);
    try std.testing.expectEqual(405, y5[0].im);
    try std.testing.expectEqual(367, y5[2].re);
    try std.testing.expectEqual(465, y5[2].im);
    try std.testing.expectEqual(298, y5[4].re);
    try std.testing.expectEqual(526, y5[4].im);
    try std.testing.expectEqual(42, y5[6].re);
    try std.testing.expectEqual(606, y5[6].im);
    try std.testing.expectEqual(-535, y5[8].re);
    try std.testing.expectEqual(735, y5[8].im);

    const x6 = try a.alloc(cf64, 2 * n);
    defer a.free(x6);
    const y6 = try a.alloc(cf64, 2 * n);
    defer a.free(y6);

    @memcpy(x6.ptr, &[_]cf64{
        cf64.init(9, 10),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y6.ptr, &[_]cf64{
        cf64.init(5, 5),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    hpmv(cf64, .RowMajor, .Lower, n, alpha, A.ptr, x6.ptr, -2, beta, y6.ptr, -2);

    try std.testing.expectEqual(347, y6[8].re);
    try std.testing.expectEqual(405, y6[8].im);
    try std.testing.expectEqual(367, y6[6].re);
    try std.testing.expectEqual(465, y6[6].im);
    try std.testing.expectEqual(298, y6[4].re);
    try std.testing.expectEqual(526, y6[4].im);
    try std.testing.expectEqual(42, y6[2].re);
    try std.testing.expectEqual(606, y6[2].im);
    try std.testing.expectEqual(-535, y6[0].re);
    try std.testing.expectEqual(735, y6[0].im);

    const x7 = try a.alloc(cf64, 2 * n);
    defer a.free(x7);
    const y7 = try a.alloc(cf64, 2 * n);
    defer a.free(y7);

    @memcpy(x7.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(9, 10),
        cf64.init(0, 0),
    });
    @memcpy(y7.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(5, 5),
        cf64.init(0, 0),
    });

    hpmv(cf64, .ColumnMajor, .Lower, n, alpha, A.ptr, x7.ptr, 2, beta, y7.ptr, 2);

    try std.testing.expectEqual(187, y7[0].re);
    try std.testing.expectEqual(225, y7[0].im);
    try std.testing.expectEqual(330, y7[2].re);
    try std.testing.expectEqual(450, y7[2].im);
    try std.testing.expectEqual(292, y7[4].re);
    try std.testing.expectEqual(592, y7[4].im);
    try std.testing.expectEqual(27, y7[6].re);
    try std.testing.expectEqual(665, y7[6].im);
    try std.testing.expectEqual(-475, y7[8].re);
    try std.testing.expectEqual(695, y7[8].im);

    const x8 = try a.alloc(cf64, 2 * n);
    defer a.free(x8);
    const y8 = try a.alloc(cf64, 2 * n);
    defer a.free(y8);

    @memcpy(x8.ptr, &[_]cf64{
        cf64.init(9, 10),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y8.ptr, &[_]cf64{
        cf64.init(5, 5),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    hpmv(cf64, .ColumnMajor, .Lower, n, alpha, A.ptr, x8.ptr, -2, beta, y8.ptr, -2);

    try std.testing.expectEqual(187, y8[8].re);
    try std.testing.expectEqual(225, y8[8].im);
    try std.testing.expectEqual(330, y8[6].re);
    try std.testing.expectEqual(450, y8[6].im);
    try std.testing.expectEqual(292, y8[4].re);
    try std.testing.expectEqual(592, y8[4].im);
    try std.testing.expectEqual(27, y8[2].re);
    try std.testing.expectEqual(665, y8[2].im);
    try std.testing.expectEqual(-475, y8[0].re);
    try std.testing.expectEqual(695, y8[0].im);
}
