const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const her = zml.linalg.blas.her;

test her {
    const a = std.testing.allocator;

    const n = 5;
    const alpha: f64 = 2;

    const A = try a.alloc(cf64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(cf64, 2 * n);
    defer a.free(x1);

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

    her(.row_major, .upper, n, alpha, x1.ptr, 2, A.ptr, n, .{}) catch unreachable;

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

    const x2 = try a.alloc(cf64, 2 * n);
    defer a.free(x2);

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

    her(.row_major, .upper, n, alpha, x2.ptr, -2, A.ptr, n, .{}) catch unreachable;

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

    const x3 = try a.alloc(cf64, 2 * n);
    defer a.free(x3);

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

    her(.col_major, .upper, n, alpha, x3.ptr, 2, A.ptr, n, .{}) catch unreachable;

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

    const x4 = try a.alloc(cf64, 2 * n);
    defer a.free(x4);

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

    her(.col_major, .upper, n, alpha, x4.ptr, -2, A.ptr, n, .{}) catch unreachable;

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

    const x5 = try a.alloc(cf64, 2 * n);
    defer a.free(x5);

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

    her(.row_major, .lower, n, alpha, x5.ptr, 2, A.ptr, n, .{}) catch unreachable;

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

    const x6 = try a.alloc(cf64, 2 * n);
    defer a.free(x6);

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

    her(.row_major, .lower, n, alpha, x6.ptr, -2, A.ptr, n, .{}) catch unreachable;

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

    const x7 = try a.alloc(cf64, 2 * n);
    defer a.free(x7);

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

    her(.col_major, .lower, n, alpha, x7.ptr, 2, A.ptr, n, .{}) catch unreachable;

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

    const x8 = try a.alloc(cf64, 2 * n);
    defer a.free(x8);

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

    her(.col_major, .lower, n, alpha, x8.ptr, -2, A.ptr, n, .{}) catch unreachable;

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
