const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const hbmv = zml.linalg.blas.hbmv;

test hbmv {
    const a = std.testing.allocator;

    const n = 5;
    const k = 1;
    const alpha = cf64.init(1, 1);
    const beta = cf64.init(3, 3);

    const A = try a.alloc(cf64, (1 + 2 * k) * n);
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

    hbmv(.row_major, .upper, n, k, alpha, A.ptr, k + 1, x1.ptr, 2, beta, y1.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(-17, y1[0].re);
    try std.testing.expectEqual(21, y1[0].im);
    try std.testing.expectEqual(-47, y1[2].re);
    try std.testing.expectEqual(81, y1[2].im);
    try std.testing.expectEqual(-77, y1[4].re);
    try std.testing.expectEqual(189, y1[4].im);
    try std.testing.expectEqual(-107, y1[6].re);
    try std.testing.expectEqual(345, y1[6].im);
    try std.testing.expectEqual(103, y1[8].re);
    try std.testing.expectEqual(329, y1[8].im);

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

    hbmv(.row_major, .upper, n, k, alpha, A.ptr, k + 1, x2.ptr, -2, beta, y2.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(-17, y2[8].re);
    try std.testing.expectEqual(21, y2[8].im);
    try std.testing.expectEqual(-47, y2[6].re);
    try std.testing.expectEqual(81, y2[6].im);
    try std.testing.expectEqual(-77, y2[4].re);
    try std.testing.expectEqual(189, y2[4].im);
    try std.testing.expectEqual(-107, y2[2].re);
    try std.testing.expectEqual(345, y2[2].im);
    try std.testing.expectEqual(103, y2[0].re);
    try std.testing.expectEqual(329, y2[0].im);

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

    hbmv(.col_major, .upper, n, k, alpha, A.ptr, k + 1, x3.ptr, 2, beta, y3.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(-26, y3[0].re);
    try std.testing.expectEqual(30, y3[0].im);
    try std.testing.expectEqual(-58, y3[2].re);
    try std.testing.expectEqual(102, y3[2].im);
    try std.testing.expectEqual(-88, y3[4].re);
    try std.testing.expectEqual(222, y3[4].im);
    try std.testing.expectEqual(-118, y3[6].re);
    try std.testing.expectEqual(390, y3[6].im);
    try std.testing.expectEqual(116, y3[8].re);
    try std.testing.expectEqual(364, y3[8].im);

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

    hbmv(.col_major, .upper, n, k, alpha, A.ptr, k + 1, x4.ptr, -2, beta, y4.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(-26, y4[8].re);
    try std.testing.expectEqual(30, y4[8].im);
    try std.testing.expectEqual(-58, y4[6].re);
    try std.testing.expectEqual(102, y4[6].im);
    try std.testing.expectEqual(-88, y4[4].re);
    try std.testing.expectEqual(222, y4[4].im);
    try std.testing.expectEqual(-118, y4[2].re);
    try std.testing.expectEqual(390, y4[2].im);
    try std.testing.expectEqual(116, y4[0].re);
    try std.testing.expectEqual(364, y4[0].im);

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

    hbmv(.row_major, .lower, n, k, alpha, A.ptr, k + 1, x5.ptr, 2, beta, y5.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(16, y5[0].re);
    try std.testing.expectEqual(36, y5[0].im);
    try std.testing.expectEqual(34, y5[2].re);
    try std.testing.expectEqual(106, y5[2].im);
    try std.testing.expectEqual(52, y5[4].re);
    try std.testing.expectEqual(226, y5[4].im);
    try std.testing.expectEqual(70, y5[6].re);
    try std.testing.expectEqual(394, y5[6].im);
    try std.testing.expectEqual(-154, y5[8].re);
    try std.testing.expectEqual(346, y5[8].im);

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

    hbmv(.row_major, .lower, n, k, alpha, A.ptr, k + 1, x6.ptr, -2, beta, y6.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(16, y6[8].re);
    try std.testing.expectEqual(36, y6[8].im);
    try std.testing.expectEqual(34, y6[6].re);
    try std.testing.expectEqual(106, y6[6].im);
    try std.testing.expectEqual(52, y6[4].re);
    try std.testing.expectEqual(226, y6[4].im);
    try std.testing.expectEqual(70, y6[2].re);
    try std.testing.expectEqual(394, y6[2].im);
    try std.testing.expectEqual(-154, y6[0].re);
    try std.testing.expectEqual(346, y6[0].im);

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

    hbmv(.col_major, .lower, n, k, alpha, A.ptr, k + 1, x7.ptr, 2, beta, y7.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(11, y7[0].re);
    try std.testing.expectEqual(25, y7[0].im);
    try std.testing.expectEqual(29, y7[2].re);
    try std.testing.expectEqual(85, y7[2].im);
    try std.testing.expectEqual(47, y7[4].re);
    try std.testing.expectEqual(193, y7[4].im);
    try std.testing.expectEqual(65, y7[6].re);
    try std.testing.expectEqual(349, y7[6].im);
    try std.testing.expectEqual(-137, y7[8].re);
    try std.testing.expectEqual(313, y7[8].im);

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

    hbmv(.col_major, .lower, n, k, alpha, A.ptr, k + 1, x8.ptr, -2, beta, y8.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(11, y8[8].re);
    try std.testing.expectEqual(25, y8[8].im);
    try std.testing.expectEqual(29, y8[6].re);
    try std.testing.expectEqual(85, y8[6].im);
    try std.testing.expectEqual(47, y8[4].re);
    try std.testing.expectEqual(193, y8[4].im);
    try std.testing.expectEqual(65, y8[2].re);
    try std.testing.expectEqual(349, y8[2].im);
    try std.testing.expectEqual(-137, y8[0].re);
    try std.testing.expectEqual(313, y8[0].im);
}
