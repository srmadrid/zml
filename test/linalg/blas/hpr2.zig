const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const hpr2 = zml.linalg.blas.hpr2;

test hpr2 {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = cf64.init(1, 1);

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

    hpr2(.row_major, .upper, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(5, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(12, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(19, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(26, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(33, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(30, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(45, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(60, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(75, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(70, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(93, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(116, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(125, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(156, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(195, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.row_major, .upper, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(9, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(22, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(35, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(48, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(61, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(54, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(83, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(112, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(141, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(130, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(175, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(220, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(237, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(298, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(375, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.col_major, .upper, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(13, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(32, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(59, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(64, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(99, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(114, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(105, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(164, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(223, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(242, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(203, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(286, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(341, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(440, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(555, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.col_major, .upper, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(17, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(42, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(83, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(80, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(137, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(174, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(127, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(216, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(305, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(354, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(231, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(352, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(445, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(582, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(735, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.row_major, .lower, n, alpha, x5.ptr, 2, y5.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(21, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(52, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(107, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(96, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(175, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(234, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(149, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(268, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(387, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(466, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(259, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(418, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(549, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(724, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(915, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.row_major, .lower, n, alpha, x6.ptr, -2, y6.ptr, -2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(25, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(62, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(131, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(112, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(213, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(294, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(171, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(320, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(469, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(578, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(287, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(484, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(653, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(866, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(1095, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.col_major, .lower, n, alpha, x7.ptr, 2, y7.ptr, 2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(29, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(72, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(147, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(134, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(241, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(318, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(209, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(372, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(535, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(638, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(369, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(588, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(765, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(1008, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(1275, A[14].re);
    try std.testing.expectEqual(0, A[14].im);

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

    hpr2(.col_major, .lower, n, alpha, x8.ptr, -2, y8.ptr, -2, A.ptr, .{}) catch unreachable;

    try std.testing.expectEqual(33, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(82, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(163, A[2].re);
    try std.testing.expectEqual(0, A[2].im);
    try std.testing.expectEqual(156, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(269, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(342, A[5].re);
    try std.testing.expectEqual(0, A[5].im);
    try std.testing.expectEqual(247, A[6].re);
    try std.testing.expectEqual(7, A[6].im);
    try std.testing.expectEqual(424, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(601, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(698, A[9].re);
    try std.testing.expectEqual(0, A[9].im);
    try std.testing.expectEqual(451, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(692, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(877, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(1150, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(1455, A[14].re);
    try std.testing.expectEqual(0, A[14].im);
}
