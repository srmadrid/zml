const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const hpr = zml.linalg.blas.hpr;

test hpr {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = 2;

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

    hpr(.row_major, .upper, n, alpha, x1.ptr, 2, A.ptr, .{}) catch unreachable;

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
    try std.testing.expectEqual(56, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(85, A[6].re);
    try std.testing.expectEqual(11, A[6].im);
    try std.testing.expectEqual(114, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(143, A[8].re);
    try std.testing.expectEqual(21, A[8].im);
    try std.testing.expectEqual(132, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(177, A[10].re);
    try std.testing.expectEqual(15, A[10].im);
    try std.testing.expectEqual(222, A[11].re);
    try std.testing.expectEqual(20, A[11].im);
    try std.testing.expectEqual(239, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(300, A[13].re);
    try std.testing.expectEqual(18, A[13].im);
    try std.testing.expectEqual(377, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.row_major, .upper, n, alpha, x2.ptr, -2, A.ptr, .{}) catch unreachable;

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
    try std.testing.expectEqual(106, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(163, A[6].re);
    try std.testing.expectEqual(15, A[6].im);
    try std.testing.expectEqual(220, A[7].re);
    try std.testing.expectEqual(24, A[7].im);
    try std.testing.expectEqual(277, A[8].re);
    try std.testing.expectEqual(33, A[8].im);
    try std.testing.expectEqual(254, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(343, A[10].re);
    try std.testing.expectEqual(19, A[10].im);
    try std.testing.expectEqual(432, A[11].re);
    try std.testing.expectEqual(28, A[11].im);
    try std.testing.expectEqual(465, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(586, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(739, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.col_major, .upper, n, alpha, x3.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(31, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(68, A[1].re);
    try std.testing.expectEqual(14, A[1].im);
    try std.testing.expectEqual(121, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(130, A[3].re);
    try std.testing.expectEqual(36, A[3].im);
    try std.testing.expectEqual(199, A[4].re);
    try std.testing.expectEqual(41, A[4].im);
    try std.testing.expectEqual(228, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(209, A[6].re);
    try std.testing.expectEqual(27, A[6].im);
    try std.testing.expectEqual(326, A[7].re);
    try std.testing.expectEqual(32, A[7].im);
    try std.testing.expectEqual(443, A[8].re);
    try std.testing.expectEqual(37, A[8].im);
    try std.testing.expectEqual(480, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(401, A[10].re);
    try std.testing.expectEqual(35, A[10].im);
    try std.testing.expectEqual(566, A[11].re);
    try std.testing.expectEqual(40, A[11].im);
    try std.testing.expectEqual(675, A[12].re);
    try std.testing.expectEqual(8, A[12].im);
    try std.testing.expectEqual(872, A[13].re);
    try std.testing.expectEqual(26, A[13].im);
    try std.testing.expectEqual(1101, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.col_major, .upper, n, alpha, x4.ptr, -2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(41, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(90, A[1].re);
    try std.testing.expectEqual(18, A[1].im);
    try std.testing.expectEqual(171, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(164, A[3].re);
    try std.testing.expectEqual(44, A[3].im);
    try std.testing.expectEqual(277, A[4].re);
    try std.testing.expectEqual(45, A[4].im);
    try std.testing.expectEqual(350, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(255, A[6].re);
    try std.testing.expectEqual(39, A[6].im);
    try std.testing.expectEqual(432, A[7].re);
    try std.testing.expectEqual(40, A[7].im);
    try std.testing.expectEqual(609, A[8].re);
    try std.testing.expectEqual(41, A[8].im);
    try std.testing.expectEqual(706, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(459, A[10].re);
    try std.testing.expectEqual(51, A[10].im);
    try std.testing.expectEqual(700, A[11].re);
    try std.testing.expectEqual(52, A[11].im);
    try std.testing.expectEqual(885, A[12].re);
    try std.testing.expectEqual(16, A[12].im);
    try std.testing.expectEqual(1158, A[13].re);
    try std.testing.expectEqual(30, A[13].im);
    try std.testing.expectEqual(1463, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.row_major, .lower, n, alpha, x5.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(51, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(112, A[1].re);
    try std.testing.expectEqual(14, A[1].im);
    try std.testing.expectEqual(221, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(198, A[3].re);
    try std.testing.expectEqual(36, A[3].im);
    try std.testing.expectEqual(355, A[4].re);
    try std.testing.expectEqual(41, A[4].im);
    try std.testing.expectEqual(472, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(301, A[6].re);
    try std.testing.expectEqual(27, A[6].im);
    try std.testing.expectEqual(538, A[7].re);
    try std.testing.expectEqual(32, A[7].im);
    try std.testing.expectEqual(775, A[8].re);
    try std.testing.expectEqual(37, A[8].im);
    try std.testing.expectEqual(932, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(517, A[10].re);
    try std.testing.expectEqual(35, A[10].im);
    try std.testing.expectEqual(834, A[11].re);
    try std.testing.expectEqual(40, A[11].im);
    try std.testing.expectEqual(1095, A[12].re);
    try std.testing.expectEqual(8, A[12].im);
    try std.testing.expectEqual(1444, A[13].re);
    try std.testing.expectEqual(26, A[13].im);
    try std.testing.expectEqual(1825, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.row_major, .lower, n, alpha, x6.ptr, -2, A.ptr, .{}) catch unreachable;
    try std.testing.expectEqual(61, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(134, A[1].re);
    try std.testing.expectEqual(10, A[1].im);
    try std.testing.expectEqual(271, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(232, A[3].re);
    try std.testing.expectEqual(28, A[3].im);
    try std.testing.expectEqual(433, A[4].re);
    try std.testing.expectEqual(37, A[4].im);
    try std.testing.expectEqual(594, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(347, A[6].re);
    try std.testing.expectEqual(15, A[6].im);
    try std.testing.expectEqual(644, A[7].re);
    try std.testing.expectEqual(24, A[7].im);
    try std.testing.expectEqual(941, A[8].re);
    try std.testing.expectEqual(33, A[8].im);
    try std.testing.expectEqual(1158, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(575, A[10].re);
    try std.testing.expectEqual(19, A[10].im);
    try std.testing.expectEqual(968, A[11].re);
    try std.testing.expectEqual(28, A[11].im);
    try std.testing.expectEqual(1305, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(1730, A[13].re);
    try std.testing.expectEqual(22, A[13].im);
    try std.testing.expectEqual(2187, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.col_major, .lower, n, alpha, x7.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(71, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(156, A[1].re);
    try std.testing.expectEqual(6, A[1].im);
    try std.testing.expectEqual(305, A[2].re);
    try std.testing.expectEqual(-8, A[2].im);
    try std.testing.expectEqual(278, A[3].re);
    try std.testing.expectEqual(16, A[3].im);
    try std.testing.expectEqual(491, A[4].re);
    try std.testing.expectEqual(21, A[4].im);
    try std.testing.expectEqual(644, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(425, A[6].re);
    try std.testing.expectEqual(11, A[6].im);
    try std.testing.expectEqual(750, A[7].re);
    try std.testing.expectEqual(16, A[7].im);
    try std.testing.expectEqual(1075, A[8].re);
    try std.testing.expectEqual(21, A[8].im);
    try std.testing.expectEqual(1280, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(741, A[10].re);
    try std.testing.expectEqual(15, A[10].im);
    try std.testing.expectEqual(1178, A[11].re);
    try std.testing.expectEqual(20, A[11].im);
    try std.testing.expectEqual(1531, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(2016, A[13].re);
    try std.testing.expectEqual(18, A[13].im);
    try std.testing.expectEqual(2549, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr(.col_major, .lower, n, alpha, x8.ptr, -2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(81, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(178, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(339, A[2].re);
    try std.testing.expectEqual(-16, A[2].im);
    try std.testing.expectEqual(324, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(549, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(694, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(503, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(856, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(1209, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(1402, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(907, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(1388, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(1757, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(2302, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(2911, A[14].re);
    try std.testing.expectEqual(0, A[14].im);
}
