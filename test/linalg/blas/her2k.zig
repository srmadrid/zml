const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const her2k = zml.linalg.blas.her2k;

test her2k {
    const a = std.testing.allocator;

    const n = 5;
    const k = 3;
    const alpha = cf64.init(2, 2);
    const beta = 3;

    const A = try a.alloc(cf64, n * k);
    defer a.free(A);
    const B = try a.alloc(cf64, n * k);
    defer a.free(B);
    const C = try a.alloc(cf64, n * n);
    defer a.free(C);

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
    @memcpy(B.ptr, &[_]cf64{
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
    @memcpy(C.ptr, &[_]cf64{
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

    her2k(cf64, .RowMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(115, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(262, C[1].re);
    try std.testing.expectEqual(6, C[1].im);
    try std.testing.expectEqual(409, C[2].re);
    try std.testing.expectEqual(9, C[2].im);
    try std.testing.expectEqual(556, C[3].re);
    try std.testing.expectEqual(12, C[3].im);
    try std.testing.expectEqual(703, C[4].re);
    try std.testing.expectEqual(15, C[4].im);
    try std.testing.expectEqual(6, C[5].re);
    try std.testing.expectEqual(6, C[5].im);
    try std.testing.expectEqual(637, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(1000, C[7].re);
    try std.testing.expectEqual(24, C[7].im);
    try std.testing.expectEqual(1363, C[8].re);
    try std.testing.expectEqual(27, C[8].im);
    try std.testing.expectEqual(1726, C[9].re);
    try std.testing.expectEqual(30, C[9].im);
    try std.testing.expectEqual(11, C[10].re);
    try std.testing.expectEqual(11, C[10].im);
    try std.testing.expectEqual(12, C[11].re);
    try std.testing.expectEqual(12, C[11].im);
    try std.testing.expectEqual(1591, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(2170, C[13].re);
    try std.testing.expectEqual(42, C[13].im);
    try std.testing.expectEqual(2749, C[14].re);
    try std.testing.expectEqual(45, C[14].im);
    try std.testing.expectEqual(16, C[15].re);
    try std.testing.expectEqual(16, C[15].im);
    try std.testing.expectEqual(17, C[16].re);
    try std.testing.expectEqual(17, C[16].im);
    try std.testing.expectEqual(18, C[17].re);
    try std.testing.expectEqual(18, C[17].im);
    try std.testing.expectEqual(2977, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(3772, C[19].re);
    try std.testing.expectEqual(60, C[19].im);
    try std.testing.expectEqual(21, C[20].re);
    try std.testing.expectEqual(21, C[20].im);
    try std.testing.expectEqual(22, C[21].re);
    try std.testing.expectEqual(22, C[21].im);
    try std.testing.expectEqual(23, C[22].re);
    try std.testing.expectEqual(23, C[22].im);
    try std.testing.expectEqual(24, C[23].re);
    try std.testing.expectEqual(24, C[23].im);
    try std.testing.expectEqual(4795, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .ColumnMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(1609, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(262, C[1].re);
    try std.testing.expectEqual(6, C[1].im);
    try std.testing.expectEqual(409, C[2].re);
    try std.testing.expectEqual(9, C[2].im);
    try std.testing.expectEqual(556, C[3].re);
    try std.testing.expectEqual(12, C[3].im);
    try std.testing.expectEqual(703, C[4].re);
    try std.testing.expectEqual(15, C[4].im);
    try std.testing.expectEqual(1426, C[5].re);
    try std.testing.expectEqual(18, C[5].im);
    try std.testing.expectEqual(3487, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(1000, C[7].re);
    try std.testing.expectEqual(24, C[7].im);
    try std.testing.expectEqual(1363, C[8].re);
    try std.testing.expectEqual(27, C[8].im);
    try std.testing.expectEqual(1726, C[9].re);
    try std.testing.expectEqual(30, C[9].im);
    try std.testing.expectEqual(1585, C[10].re);
    try std.testing.expectEqual(33, C[10].im);
    try std.testing.expectEqual(1780, C[11].re);
    try std.testing.expectEqual(36, C[11].im);
    try std.testing.expectEqual(6709, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(2170, C[13].re);
    try std.testing.expectEqual(42, C[13].im);
    try std.testing.expectEqual(2749, C[14].re);
    try std.testing.expectEqual(45, C[14].im);
    try std.testing.expectEqual(1744, C[15].re);
    try std.testing.expectEqual(48, C[15].im);
    try std.testing.expectEqual(1963, C[16].re);
    try std.testing.expectEqual(51, C[16].im);
    try std.testing.expectEqual(2182, C[17].re);
    try std.testing.expectEqual(54, C[17].im);
    try std.testing.expectEqual(11275, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(3772, C[19].re);
    try std.testing.expectEqual(60, C[19].im);
    try std.testing.expectEqual(1903, C[20].re);
    try std.testing.expectEqual(63, C[20].im);
    try std.testing.expectEqual(2146, C[21].re);
    try std.testing.expectEqual(66, C[21].im);
    try std.testing.expectEqual(2389, C[22].re);
    try std.testing.expectEqual(69, C[22].im);
    try std.testing.expectEqual(2632, C[23].re);
    try std.testing.expectEqual(72, C[23].im);
    try std.testing.expectEqual(17185, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .RowMajor, .Upper, .ConjTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(6091, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(2194, C[1].re);
    try std.testing.expectEqual(18, C[1].im);
    try std.testing.expectEqual(2779, C[2].re);
    try std.testing.expectEqual(27, C[2].im);
    try std.testing.expectEqual(3364, C[3].re);
    try std.testing.expectEqual(36, C[3].im);
    try std.testing.expectEqual(3949, C[4].re);
    try std.testing.expectEqual(45, C[4].im);
    try std.testing.expectEqual(1426, C[5].re);
    try std.testing.expectEqual(18, C[5].im);
    try std.testing.expectEqual(12037, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(4744, C[7].re);
    try std.testing.expectEqual(72, C[7].im);
    try std.testing.expectEqual(6001, C[8].re);
    try std.testing.expectEqual(81, C[8].im);
    try std.testing.expectEqual(7258, C[9].re);
    try std.testing.expectEqual(90, C[9].im);
    try std.testing.expectEqual(1585, C[10].re);
    try std.testing.expectEqual(33, C[10].im);
    try std.testing.expectEqual(1780, C[11].re);
    try std.testing.expectEqual(36, C[11].im);
    try std.testing.expectEqual(22063, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(8638, C[13].re);
    try std.testing.expectEqual(126, C[13].im);
    try std.testing.expectEqual(10567, C[14].re);
    try std.testing.expectEqual(135, C[14].im);
    try std.testing.expectEqual(1744, C[15].re);
    try std.testing.expectEqual(48, C[15].im);
    try std.testing.expectEqual(1963, C[16].re);
    try std.testing.expectEqual(51, C[16].im);
    try std.testing.expectEqual(2182, C[17].re);
    try std.testing.expectEqual(54, C[17].im);
    try std.testing.expectEqual(36169, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(13876, C[19].re);
    try std.testing.expectEqual(180, C[19].im);
    try std.testing.expectEqual(1903, C[20].re);
    try std.testing.expectEqual(63, C[20].im);
    try std.testing.expectEqual(2146, C[21].re);
    try std.testing.expectEqual(66, C[21].im);
    try std.testing.expectEqual(2389, C[22].re);
    try std.testing.expectEqual(69, C[22].im);
    try std.testing.expectEqual(2632, C[23].re);
    try std.testing.expectEqual(72, C[23].im);
    try std.testing.expectEqual(54355, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .ColumnMajor, .Upper, .ConjTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(18385, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(2194, C[1].re);
    try std.testing.expectEqual(18, C[1].im);
    try std.testing.expectEqual(2779, C[2].re);
    try std.testing.expectEqual(27, C[2].im);
    try std.testing.expectEqual(3364, C[3].re);
    try std.testing.expectEqual(36, C[3].im);
    try std.testing.expectEqual(3949, C[4].re);
    try std.testing.expectEqual(45, C[4].im);
    try std.testing.expectEqual(4534, C[5].re);
    try std.testing.expectEqual(54, C[5].im);
    try std.testing.expectEqual(36727, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(4744, C[7].re);
    try std.testing.expectEqual(72, C[7].im);
    try std.testing.expectEqual(6001, C[8].re);
    try std.testing.expectEqual(81, C[8].im);
    try std.testing.expectEqual(7258, C[9].re);
    try std.testing.expectEqual(90, C[9].im);
    try std.testing.expectEqual(5155, C[10].re);
    try std.testing.expectEqual(99, C[10].im);
    try std.testing.expectEqual(6316, C[11].re);
    try std.testing.expectEqual(108, C[11].im);
    try std.testing.expectEqual(67741, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(8638, C[13].re);
    try std.testing.expectEqual(126, C[13].im);
    try std.testing.expectEqual(10567, C[14].re);
    try std.testing.expectEqual(135, C[14].im);
    try std.testing.expectEqual(5776, C[15].re);
    try std.testing.expectEqual(144, C[15].im);
    try std.testing.expectEqual(7225, C[16].re);
    try std.testing.expectEqual(153, C[16].im);
    try std.testing.expectEqual(8674, C[17].re);
    try std.testing.expectEqual(162, C[17].im);
    try std.testing.expectEqual(111427, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(13876, C[19].re);
    try std.testing.expectEqual(180, C[19].im);
    try std.testing.expectEqual(6397, C[20].re);
    try std.testing.expectEqual(189, C[20].im);
    try std.testing.expectEqual(8134, C[21].re);
    try std.testing.expectEqual(198, C[21].im);
    try std.testing.expectEqual(9871, C[22].re);
    try std.testing.expectEqual(207, C[22].im);
    try std.testing.expectEqual(11608, C[23].re);
    try std.testing.expectEqual(216, C[23].im);
    try std.testing.expectEqual(167785, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .RowMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(55267, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(2194, C[1].re);
    try std.testing.expectEqual(18, C[1].im);
    try std.testing.expectEqual(2779, C[2].re);
    try std.testing.expectEqual(27, C[2].im);
    try std.testing.expectEqual(3364, C[3].re);
    try std.testing.expectEqual(36, C[3].im);
    try std.testing.expectEqual(3949, C[4].re);
    try std.testing.expectEqual(45, C[4].im);
    try std.testing.expectEqual(13858, C[5].re);
    try std.testing.expectEqual(162, C[5].im);
    try std.testing.expectEqual(110797, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(4744, C[7].re);
    try std.testing.expectEqual(72, C[7].im);
    try std.testing.expectEqual(6001, C[8].re);
    try std.testing.expectEqual(81, C[8].im);
    try std.testing.expectEqual(7258, C[9].re);
    try std.testing.expectEqual(90, C[9].im);
    try std.testing.expectEqual(15865, C[10].re);
    try std.testing.expectEqual(297, C[10].im);
    try std.testing.expectEqual(19924, C[11].re);
    try std.testing.expectEqual(324, C[11].im);
    try std.testing.expectEqual(204775, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(8638, C[13].re);
    try std.testing.expectEqual(126, C[13].im);
    try std.testing.expectEqual(10567, C[14].re);
    try std.testing.expectEqual(135, C[14].im);
    try std.testing.expectEqual(17872, C[15].re);
    try std.testing.expectEqual(432, C[15].im);
    try std.testing.expectEqual(23011, C[16].re);
    try std.testing.expectEqual(459, C[16].im);
    try std.testing.expectEqual(28150, C[17].re);
    try std.testing.expectEqual(486, C[17].im);
    try std.testing.expectEqual(337201, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(13876, C[19].re);
    try std.testing.expectEqual(180, C[19].im);
    try std.testing.expectEqual(19879, C[20].re);
    try std.testing.expectEqual(567, C[20].im);
    try std.testing.expectEqual(26098, C[21].re);
    try std.testing.expectEqual(594, C[21].im);
    try std.testing.expectEqual(32317, C[22].re);
    try std.testing.expectEqual(621, C[22].im);
    try std.testing.expectEqual(38536, C[23].re);
    try std.testing.expectEqual(648, C[23].im);
    try std.testing.expectEqual(508075, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .ColumnMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(167065, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(7990, C[1].re);
    try std.testing.expectEqual(54, C[1].im);
    try std.testing.expectEqual(9889, C[2].re);
    try std.testing.expectEqual(81, C[2].im);
    try std.testing.expectEqual(11788, C[3].re);
    try std.testing.expectEqual(108, C[3].im);
    try std.testing.expectEqual(13687, C[4].re);
    try std.testing.expectEqual(135, C[4].im);
    try std.testing.expectEqual(13858, C[5].re);
    try std.testing.expectEqual(162, C[5].im);
    try std.testing.expectEqual(333967, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(15976, C[7].re);
    try std.testing.expectEqual(216, C[7].im);
    try std.testing.expectEqual(19915, C[8].re);
    try std.testing.expectEqual(243, C[8].im);
    try std.testing.expectEqual(23854, C[9].re);
    try std.testing.expectEqual(270, C[9].im);
    try std.testing.expectEqual(15865, C[10].re);
    try std.testing.expectEqual(297, C[10].im);
    try std.testing.expectEqual(19924, C[11].re);
    try std.testing.expectEqual(324, C[11].im);
    try std.testing.expectEqual(616261, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(28042, C[13].re);
    try std.testing.expectEqual(378, C[13].im);
    try std.testing.expectEqual(34021, C[14].re);
    try std.testing.expectEqual(405, C[14].im);
    try std.testing.expectEqual(17872, C[15].re);
    try std.testing.expectEqual(432, C[15].im);
    try std.testing.expectEqual(23011, C[16].re);
    try std.testing.expectEqual(459, C[16].im);
    try std.testing.expectEqual(28150, C[17].re);
    try std.testing.expectEqual(486, C[17].im);
    try std.testing.expectEqual(1013947, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(44188, C[19].re);
    try std.testing.expectEqual(540, C[19].im);
    try std.testing.expectEqual(19879, C[20].re);
    try std.testing.expectEqual(567, C[20].im);
    try std.testing.expectEqual(26098, C[21].re);
    try std.testing.expectEqual(594, C[21].im);
    try std.testing.expectEqual(32317, C[22].re);
    try std.testing.expectEqual(621, C[22].im);
    try std.testing.expectEqual(38536, C[23].re);
    try std.testing.expectEqual(648, C[23].im);
    try std.testing.expectEqual(1527025, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .RowMajor, .Lower, .ConjTrans, n, k, alpha, A.ptr, n, B.ptr, n, beta, C.ptr, n);

    try std.testing.expectEqual(502459, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(7990, C[1].re);
    try std.testing.expectEqual(54, C[1].im);
    try std.testing.expectEqual(9889, C[2].re);
    try std.testing.expectEqual(81, C[2].im);
    try std.testing.expectEqual(11788, C[3].re);
    try std.testing.expectEqual(108, C[3].im);
    try std.testing.expectEqual(13687, C[4].re);
    try std.testing.expectEqual(135, C[4].im);
    try std.testing.expectEqual(42982, C[5].re);
    try std.testing.expectEqual(486, C[5].im);
    try std.testing.expectEqual(1003477, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(15976, C[7].re);
    try std.testing.expectEqual(216, C[7].im);
    try std.testing.expectEqual(19915, C[8].re);
    try std.testing.expectEqual(243, C[8].im);
    try std.testing.expectEqual(23854, C[9].re);
    try std.testing.expectEqual(270, C[9].im);
    try std.testing.expectEqual(49147, C[10].re);
    try std.testing.expectEqual(891, C[10].im);
    try std.testing.expectEqual(61516, C[11].re);
    try std.testing.expectEqual(972, C[11].im);
    try std.testing.expectEqual(1850719, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(28042, C[13].re);
    try std.testing.expectEqual(378, C[13].im);
    try std.testing.expectEqual(34021, C[14].re);
    try std.testing.expectEqual(405, C[14].im);
    try std.testing.expectEqual(55312, C[15].re);
    try std.testing.expectEqual(1296, C[15].im);
    try std.testing.expectEqual(70945, C[16].re);
    try std.testing.expectEqual(1377, C[16].im);
    try std.testing.expectEqual(86578, C[17].re);
    try std.testing.expectEqual(1458, C[17].im);
    try std.testing.expectEqual(3044185, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(44188, C[19].re);
    try std.testing.expectEqual(540, C[19].im);
    try std.testing.expectEqual(61477, C[20].re);
    try std.testing.expectEqual(1701, C[20].im);
    try std.testing.expectEqual(80374, C[21].re);
    try std.testing.expectEqual(1782, C[21].im);
    try std.testing.expectEqual(99271, C[22].re);
    try std.testing.expectEqual(1863, C[22].im);
    try std.testing.expectEqual(118168, C[23].re);
    try std.testing.expectEqual(1944, C[23].im);
    try std.testing.expectEqual(4583875, C[24].re);
    try std.testing.expectEqual(0, C[24].im);

    her2k(cf64, .ColumnMajor, .Lower, .ConjTrans, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n);

    try std.testing.expectEqual(1507489, C[0].re);
    try std.testing.expectEqual(0, C[0].im);
    try std.testing.expectEqual(24226, C[1].re);
    try std.testing.expectEqual(162, C[1].im);
    try std.testing.expectEqual(30067, C[2].re);
    try std.testing.expectEqual(243, C[2].im);
    try std.testing.expectEqual(35908, C[3].re);
    try std.testing.expectEqual(324, C[3].im);
    try std.testing.expectEqual(41749, C[4].re);
    try std.testing.expectEqual(405, C[4].im);
    try std.testing.expectEqual(42982, C[5].re);
    try std.testing.expectEqual(486, C[5].im);
    try std.testing.expectEqual(3011047, C[6].re);
    try std.testing.expectEqual(0, C[6].im);
    try std.testing.expectEqual(48904, C[7].re);
    try std.testing.expectEqual(648, C[7].im);
    try std.testing.expectEqual(61081, C[8].re);
    try std.testing.expectEqual(729, C[8].im);
    try std.testing.expectEqual(73258, C[9].re);
    try std.testing.expectEqual(810, C[9].im);
    try std.testing.expectEqual(49147, C[10].re);
    try std.testing.expectEqual(891, C[10].im);
    try std.testing.expectEqual(61516, C[11].re);
    try std.testing.expectEqual(972, C[11].im);
    try std.testing.expectEqual(5553709, C[12].re);
    try std.testing.expectEqual(0, C[12].im);
    try std.testing.expectEqual(86254, C[13].re);
    try std.testing.expectEqual(1134, C[13].im);
    try std.testing.expectEqual(104767, C[14].re);
    try std.testing.expectEqual(1215, C[14].im);
    try std.testing.expectEqual(55312, C[15].re);
    try std.testing.expectEqual(1296, C[15].im);
    try std.testing.expectEqual(70945, C[16].re);
    try std.testing.expectEqual(1377, C[16].im);
    try std.testing.expectEqual(86578, C[17].re);
    try std.testing.expectEqual(1458, C[17].im);
    try std.testing.expectEqual(9135475, C[18].re);
    try std.testing.expectEqual(0, C[18].im);
    try std.testing.expectEqual(136276, C[19].re);
    try std.testing.expectEqual(1620, C[19].im);
    try std.testing.expectEqual(61477, C[20].re);
    try std.testing.expectEqual(1701, C[20].im);
    try std.testing.expectEqual(80374, C[21].re);
    try std.testing.expectEqual(1782, C[21].im);
    try std.testing.expectEqual(99271, C[22].re);
    try std.testing.expectEqual(1863, C[22].im);
    try std.testing.expectEqual(118168, C[23].re);
    try std.testing.expectEqual(1944, C[23].im);
    try std.testing.expectEqual(13756345, C[24].re);
    try std.testing.expectEqual(0, C[24].im);
}
