const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const her2 = zml.linalg.blas.her2;

test her2 {
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

    her2(cf64, .RowMajor, .Upper, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr, n);

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
    try std.testing.expectEqual(6, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(31, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(46, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(61, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(76, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(11, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(12, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(73, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(96, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(119, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(16, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(17, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(18, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(131, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(162, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(21, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(22, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(23, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(24, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(205, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .RowMajor, .Upper, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr, n);

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
    try std.testing.expectEqual(6, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(55, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(84, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(113, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(142, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(11, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(12, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(133, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(178, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(223, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(16, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(17, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(18, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(243, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(304, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(21, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(22, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(23, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(24, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(385, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .ColumnMajor, .Upper, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(13, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(22, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(35, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(48, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(61, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(16, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(79, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(84, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(113, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(142, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(27, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(50, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(193, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(178, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(223, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(38, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(69, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(100, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(355, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(304, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(49, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(88, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(127, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(166, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(565, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .ColumnMajor, .Upper, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(17, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(22, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(35, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(48, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(61, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(26, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(103, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(84, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(113, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(142, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(43, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(88, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(253, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(178, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(223, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(60, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(121, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(182, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(467, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(304, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(77, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(154, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(231, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(308, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(745, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .RowMajor, .Lower, n, alpha, x5.ptr, 2, y5.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(21, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(22, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(35, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(48, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(61, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(36, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(127, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(84, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(113, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(142, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(59, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(126, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(313, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(178, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(223, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(82, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(173, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(264, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(579, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(304, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(105, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(220, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(335, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(450, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(925, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .RowMajor, .Lower, n, alpha, x6.ptr, -2, y6.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(25, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(22, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(35, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(48, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(61, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(46, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(151, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(84, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(113, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(142, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(75, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(164, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(373, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(178, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(223, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(104, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(225, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(346, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(691, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(304, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(133, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(286, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(439, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(592, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(1105, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .ColumnMajor, .Lower, n, alpha, x7.ptr, 2, y7.ptr, 2, A.ptr, n);

    try std.testing.expectEqual(29, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(32, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(51, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(70, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(89, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(46, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(175, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(122, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(165, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(208, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(75, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(164, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(433, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(260, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(327, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(104, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(225, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(346, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(803, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(446, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(133, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(286, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(439, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(592, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(1285, A[24].re);
    try std.testing.expectEqual(0, A[24].im);

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

    her2(cf64, .ColumnMajor, .Lower, n, alpha, x8.ptr, -2, y8.ptr, -2, A.ptr, n);

    try std.testing.expectEqual(33, A[0].re);
    try std.testing.expectEqual(0, A[0].im);
    try std.testing.expectEqual(42, A[1].re);
    try std.testing.expectEqual(2, A[1].im);
    try std.testing.expectEqual(67, A[2].re);
    try std.testing.expectEqual(3, A[2].im);
    try std.testing.expectEqual(92, A[3].re);
    try std.testing.expectEqual(4, A[3].im);
    try std.testing.expectEqual(117, A[4].re);
    try std.testing.expectEqual(5, A[4].im);
    try std.testing.expectEqual(46, A[5].re);
    try std.testing.expectEqual(6, A[5].im);
    try std.testing.expectEqual(199, A[6].re);
    try std.testing.expectEqual(0, A[6].im);
    try std.testing.expectEqual(160, A[7].re);
    try std.testing.expectEqual(8, A[7].im);
    try std.testing.expectEqual(217, A[8].re);
    try std.testing.expectEqual(9, A[8].im);
    try std.testing.expectEqual(274, A[9].re);
    try std.testing.expectEqual(10, A[9].im);
    try std.testing.expectEqual(75, A[10].re);
    try std.testing.expectEqual(11, A[10].im);
    try std.testing.expectEqual(164, A[11].re);
    try std.testing.expectEqual(12, A[11].im);
    try std.testing.expectEqual(493, A[12].re);
    try std.testing.expectEqual(0, A[12].im);
    try std.testing.expectEqual(342, A[13].re);
    try std.testing.expectEqual(14, A[13].im);
    try std.testing.expectEqual(431, A[14].re);
    try std.testing.expectEqual(15, A[14].im);
    try std.testing.expectEqual(104, A[15].re);
    try std.testing.expectEqual(16, A[15].im);
    try std.testing.expectEqual(225, A[16].re);
    try std.testing.expectEqual(17, A[16].im);
    try std.testing.expectEqual(346, A[17].re);
    try std.testing.expectEqual(18, A[17].im);
    try std.testing.expectEqual(915, A[18].re);
    try std.testing.expectEqual(0, A[18].im);
    try std.testing.expectEqual(588, A[19].re);
    try std.testing.expectEqual(20, A[19].im);
    try std.testing.expectEqual(133, A[20].re);
    try std.testing.expectEqual(21, A[20].im);
    try std.testing.expectEqual(286, A[21].re);
    try std.testing.expectEqual(22, A[21].im);
    try std.testing.expectEqual(439, A[22].re);
    try std.testing.expectEqual(23, A[22].im);
    try std.testing.expectEqual(592, A[23].re);
    try std.testing.expectEqual(24, A[23].im);
    try std.testing.expectEqual(1465, A[24].re);
    try std.testing.expectEqual(0, A[24].im);
}
