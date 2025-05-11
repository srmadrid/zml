const std = @import("std");
const zml = @import("zml");
const spr2 = zml.linalg.blas.spr2;

test spr2 {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = 2;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, 2 * n);
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
    @memcpy(y1.ptr, &[_]f64{
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

    spr2(f64, .RowMajor, .Upper, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr);

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

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * n);
    defer a.free(y2);

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
    @memcpy(y2.ptr, &[_]f64{
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

    spr2(f64, .RowMajor, .Upper, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr);

    try std.testing.expectEqual(9, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(45, A[4]);
    try std.testing.expectEqual(38, A[5]);
    try std.testing.expectEqual(55, A[6]);
    try std.testing.expectEqual(72, A[7]);
    try std.testing.expectEqual(89, A[8]);
    try std.testing.expectEqual(82, A[9]);
    try std.testing.expectEqual(107, A[10]);
    try std.testing.expectEqual(132, A[11]);
    try std.testing.expectEqual(141, A[12]);
    try std.testing.expectEqual(174, A[13]);
    try std.testing.expectEqual(215, A[14]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
    defer a.free(y3);

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
    @memcpy(y3.ptr, &[_]f64{
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

    spr2(f64, .ColumnMajor, .Upper, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr);

    try std.testing.expectEqual(13, A[0]);
    try std.testing.expectEqual(26, A[1]);
    try std.testing.expectEqual(43, A[2]);
    try std.testing.expectEqual(48, A[3]);
    try std.testing.expectEqual(69, A[4]);
    try std.testing.expectEqual(74, A[5]);
    try std.testing.expectEqual(71, A[6]);
    try std.testing.expectEqual(104, A[7]);
    try std.testing.expectEqual(137, A[8]);
    try std.testing.expectEqual(146, A[9]);
    try std.testing.expectEqual(127, A[10]);
    try std.testing.expectEqual(172, A[11]);
    try std.testing.expectEqual(201, A[12]);
    try std.testing.expectEqual(254, A[13]);
    try std.testing.expectEqual(315, A[14]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(f64, 2 * n);
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

    spr2(f64, .ColumnMajor, .Upper, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr);

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

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);
    const y5 = try a.alloc(f64, 2 * n);
    defer a.free(y5);

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
    @memcpy(y5.ptr, &[_]f64{
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

    spr2(f64, .RowMajor, .Lower, n, alpha, x5.ptr, 2, y5.ptr, 2, A.ptr);

    try std.testing.expectEqual(21, A[0]);
    try std.testing.expectEqual(42, A[1]);
    try std.testing.expectEqual(75, A[2]);
    try std.testing.expectEqual(72, A[3]);
    try std.testing.expectEqual(117, A[4]);
    try std.testing.expectEqual(146, A[5]);
    try std.testing.expectEqual(103, A[6]);
    try std.testing.expectEqual(168, A[7]);
    try std.testing.expectEqual(233, A[8]);
    try std.testing.expectEqual(274, A[9]);
    try std.testing.expectEqual(167, A[10]);
    try std.testing.expectEqual(252, A[11]);
    try std.testing.expectEqual(321, A[12]);
    try std.testing.expectEqual(414, A[13]);
    try std.testing.expectEqual(515, A[14]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);
    const y6 = try a.alloc(f64, 2 * n);
    defer a.free(y6);

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
    @memcpy(y6.ptr, &[_]f64{
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

    spr2(f64, .RowMajor, .Lower, n, alpha, x6.ptr, -2, y6.ptr, -2, A.ptr);

    try std.testing.expectEqual(25, A[0]);
    try std.testing.expectEqual(50, A[1]);
    try std.testing.expectEqual(91, A[2]);
    try std.testing.expectEqual(84, A[3]);
    try std.testing.expectEqual(141, A[4]);
    try std.testing.expectEqual(182, A[5]);
    try std.testing.expectEqual(119, A[6]);
    try std.testing.expectEqual(200, A[7]);
    try std.testing.expectEqual(281, A[8]);
    try std.testing.expectEqual(338, A[9]);
    try std.testing.expectEqual(187, A[10]);
    try std.testing.expectEqual(292, A[11]);
    try std.testing.expectEqual(381, A[12]);
    try std.testing.expectEqual(494, A[13]);
    try std.testing.expectEqual(615, A[14]);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);
    const y7 = try a.alloc(f64, 2 * n);
    defer a.free(y7);

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
    @memcpy(y7.ptr, &[_]f64{
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

    spr2(f64, .ColumnMajor, .Lower, n, alpha, x7.ptr, 2, y7.ptr, 2, A.ptr);

    try std.testing.expectEqual(29, A[0]);
    try std.testing.expectEqual(58, A[1]);
    try std.testing.expectEqual(103, A[2]);
    try std.testing.expectEqual(100, A[3]);
    try std.testing.expectEqual(161, A[4]);
    try std.testing.expectEqual(198, A[5]);
    try std.testing.expectEqual(143, A[6]);
    try std.testing.expectEqual(232, A[7]);
    try std.testing.expectEqual(321, A[8]);
    try std.testing.expectEqual(374, A[9]);
    try std.testing.expectEqual(235, A[10]);
    try std.testing.expectEqual(352, A[11]);
    try std.testing.expectEqual(445, A[12]);
    try std.testing.expectEqual(574, A[13]);
    try std.testing.expectEqual(715, A[14]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);
    const y8 = try a.alloc(f64, 2 * n);
    defer a.free(y8);

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
    @memcpy(y8.ptr, &[_]f64{
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

    spr2(f64, .ColumnMajor, .Lower, n, alpha, x8.ptr, -2, y8.ptr, -2, A.ptr);

    try std.testing.expectEqual(33, A[0]);
    try std.testing.expectEqual(66, A[1]);
    try std.testing.expectEqual(115, A[2]);
    try std.testing.expectEqual(116, A[3]);
    try std.testing.expectEqual(181, A[4]);
    try std.testing.expectEqual(214, A[5]);
    try std.testing.expectEqual(167, A[6]);
    try std.testing.expectEqual(264, A[7]);
    try std.testing.expectEqual(361, A[8]);
    try std.testing.expectEqual(410, A[9]);
    try std.testing.expectEqual(283, A[10]);
    try std.testing.expectEqual(412, A[11]);
    try std.testing.expectEqual(509, A[12]);
    try std.testing.expectEqual(654, A[13]);
    try std.testing.expectEqual(815, A[14]);
}
