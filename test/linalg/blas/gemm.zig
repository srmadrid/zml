const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const gemm = zml.linalg.blas.gemm;

test gemm {
    @setEvalBranchQuota(2000);
    const a = std.testing.allocator;

    const m = 4;
    const n = 5;
    const k = 3;
    const alpha: f64 = 2;
    const beta: f64 = 3;

    const A = try a.alloc(f64, m * k);
    defer a.free(A);
    const B = try a.alloc(f64, k * n);
    defer a.free(B);
    const C = try a.alloc(f64, m * n);
    defer a.free(C);

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
    });
    @memcpy(B.ptr, &[_]f64{
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
    });
    @memcpy(C.ptr, &[_]f64{
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
    });

    gemm(.row_major, .no_trans, .no_trans, m, n, k, alpha, A.ptr, k, B.ptr, n, beta, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(95, C[0]);
    try std.testing.expectEqual(110, C[1]);
    try std.testing.expectEqual(125, C[2]);
    try std.testing.expectEqual(140, C[3]);
    try std.testing.expectEqual(155, C[4]);
    try std.testing.expectEqual(218, C[5]);
    try std.testing.expectEqual(251, C[6]);
    try std.testing.expectEqual(284, C[7]);
    try std.testing.expectEqual(317, C[8]);
    try std.testing.expectEqual(350, C[9]);
    try std.testing.expectEqual(341, C[10]);
    try std.testing.expectEqual(392, C[11]);
    try std.testing.expectEqual(443, C[12]);
    try std.testing.expectEqual(494, C[13]);
    try std.testing.expectEqual(545, C[14]);
    try std.testing.expectEqual(464, C[15]);
    try std.testing.expectEqual(533, C[16]);
    try std.testing.expectEqual(602, C[17]);
    try std.testing.expectEqual(671, C[18]);
    try std.testing.expectEqual(740, C[19]);

    gemm(.col_major, .no_trans, .no_trans, m, n, k, alpha, A.ptr, m, B.ptr, k, beta, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(361, C[0]);
    try std.testing.expectEqual(418, C[1]);
    try std.testing.expectEqual(475, C[2]);
    try std.testing.expectEqual(532, C[3]);
    try std.testing.expectEqual(631, C[4]);
    try std.testing.expectEqual(850, C[5]);
    try std.testing.expectEqual(979, C[6]);
    try std.testing.expectEqual(1108, C[7]);
    try std.testing.expectEqual(1207, C[8]);
    try std.testing.expectEqual(1354, C[9]);
    try std.testing.expectEqual(1375, C[10]);
    try std.testing.expectEqual(1576, C[11]);
    try std.testing.expectEqual(1675, C[12]);
    try std.testing.expectEqual(1894, C[13]);
    try std.testing.expectEqual(2113, C[14]);
    try std.testing.expectEqual(1936, C[15]);
    try std.testing.expectEqual(2035, C[16]);
    try std.testing.expectEqual(2326, C[17]);
    try std.testing.expectEqual(2617, C[18]);
    try std.testing.expectEqual(2908, C[19]);

    gemm(.row_major, .trans, .no_trans, m, n, k, alpha, A.ptr, m, B.ptr, n, beta, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(1343, C[0]);
    try std.testing.expectEqual(1544, C[1]);
    try std.testing.expectEqual(1745, C[2]);
    try std.testing.expectEqual(1946, C[3]);
    try std.testing.expectEqual(2273, C[4]);
    try std.testing.expectEqual(2846, C[5]);
    try std.testing.expectEqual(3269, C[6]);
    try std.testing.expectEqual(3692, C[7]);
    try std.testing.expectEqual(4025, C[8]);
    try std.testing.expectEqual(4502, C[9]);
    try std.testing.expectEqual(4457, C[10]);
    try std.testing.expectEqual(5102, C[11]);
    try std.testing.expectEqual(5441, C[12]);
    try std.testing.expectEqual(6140, C[13]);
    try std.testing.expectEqual(6839, C[14]);
    try std.testing.expectEqual(6176, C[15]);
    try std.testing.expectEqual(6521, C[16]);
    try std.testing.expectEqual(7442, C[17]);
    try std.testing.expectEqual(8363, C[18]);
    try std.testing.expectEqual(9284, C[19]);

    gemm(.col_major, .trans, .no_trans, m, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(4057, C[0]);
    try std.testing.expectEqual(4696, C[1]);
    try std.testing.expectEqual(5335, C[2]);
    try std.testing.expectEqual(5974, C[3]);
    try std.testing.expectEqual(6883, C[4]);
    try std.testing.expectEqual(8692, C[5]);
    try std.testing.expectEqual(10051, C[6]);
    try std.testing.expectEqual(11410, C[7]);
    try std.testing.expectEqual(12175, C[8]);
    try std.testing.expectEqual(13750, C[9]);
    try std.testing.expectEqual(13759, C[10]);
    try std.testing.expectEqual(15838, C[11]);
    try std.testing.expectEqual(16459, C[12]);
    try std.testing.expectEqual(18754, C[13]);
    try std.testing.expectEqual(21049, C[14]);
    try std.testing.expectEqual(19258, C[15]);
    try std.testing.expectEqual(19735, C[16]);
    try std.testing.expectEqual(22750, C[17]);
    try std.testing.expectEqual(25765, C[18]);
    try std.testing.expectEqual(28780, C[19]);

    gemm(.row_major, .no_trans, .trans, m, n, k, alpha, A.ptr, k, B.ptr, k, beta, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(12199, C[0]);
    try std.testing.expectEqual(14152, C[1]);
    try std.testing.expectEqual(16105, C[2]);
    try std.testing.expectEqual(18058, C[3]);
    try std.testing.expectEqual(20821, C[4]);
    try std.testing.expectEqual(26140, C[5]);
    try std.testing.expectEqual(30307, C[6]);
    try std.testing.expectEqual(34474, C[7]);
    try std.testing.expectEqual(36859, C[8]);
    try std.testing.expectEqual(41674, C[9]);
    try std.testing.expectEqual(41377, C[10]);
    try std.testing.expectEqual(47758, C[11]);
    try std.testing.expectEqual(49765, C[12]);
    try std.testing.expectEqual(56794, C[13]);
    try std.testing.expectEqual(63823, C[14]);
    try std.testing.expectEqual(57910, C[15]);
    try std.testing.expectEqual(59539, C[16]);
    try std.testing.expectEqual(68782, C[17]);
    try std.testing.expectEqual(78025, C[18]);
    try std.testing.expectEqual(87268, C[19]);

    gemm(.col_major, .no_trans, .trans, m, n, k, alpha, A.ptr, m, B.ptr, n, beta, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(36857, C[0]);
    try std.testing.expectEqual(42752, C[1]);
    try std.testing.expectEqual(48647, C[2]);
    try std.testing.expectEqual(54542, C[3]);
    try std.testing.expectEqual(62753, C[4]);
    try std.testing.expectEqual(78752, C[5]);
    try std.testing.expectEqual(91295, C[6]);
    try std.testing.expectEqual(103838, C[7]);
    try std.testing.expectEqual(110897, C[8]);
    try std.testing.expectEqual(125390, C[9]);
    try std.testing.expectEqual(124547, C[10]);
    try std.testing.expectEqual(143738, C[11]);
    try std.testing.expectEqual(149645, C[12]);
    try std.testing.expectEqual(170786, C[13]);
    try std.testing.expectEqual(191927, C[14]);
    try std.testing.expectEqual(174242, C[15]);
    try std.testing.expectEqual(178997, C[16]);
    try std.testing.expectEqual(206786, C[17]);
    try std.testing.expectEqual(234575, C[18]);
    try std.testing.expectEqual(262364, C[19]);

    gemm(.row_major, .trans, .trans, m, n, k, alpha, A.ptr, m, B.ptr, k, beta, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(110647, C[0]);
    try std.testing.expectEqual(128422, C[1]);
    try std.testing.expectEqual(146197, C[2]);
    try std.testing.expectEqual(163972, C[3]);
    try std.testing.expectEqual(188695, C[4]);
    try std.testing.expectEqual(236344, C[5]);
    try std.testing.expectEqual(274081, C[6]);
    try std.testing.expectEqual(311818, C[7]);
    try std.testing.expectEqual(333103, C[8]);
    try std.testing.expectEqual(376690, C[9]);
    try std.testing.expectEqual(373741, C[10]);
    try std.testing.expectEqual(431440, C[11]);
    try std.testing.expectEqual(449287, C[12]);
    try std.testing.expectEqual(512836, C[13]);
    try std.testing.expectEqual(576385, C[14]);
    try std.testing.expectEqual(522838, C[15]);
    try std.testing.expectEqual(537247, C[16]);
    try std.testing.expectEqual(620758, C[17]);
    try std.testing.expectEqual(704269, C[18]);
    try std.testing.expectEqual(787780, C[19]);

    gemm(.col_major, .trans, .trans, m, n, k, alpha, A.ptr, k, B.ptr, n, beta, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(332033, C[0]);
    try std.testing.expectEqual(385466, C[1]);
    try std.testing.expectEqual(438899, C[2]);
    try std.testing.expectEqual(492332, C[3]);
    try std.testing.expectEqual(566189, C[4]);
    try std.testing.expectEqual(709262, C[5]);
    try std.testing.expectEqual(822599, C[6]);
    try std.testing.expectEqual(935936, C[7]);
    try std.testing.expectEqual(999425, C[8]);
    try std.testing.expectEqual(1130330, C[9]);
    try std.testing.expectEqual(1121627, C[10]);
    try std.testing.expectEqual(1294868, C[11]);
    try std.testing.expectEqual(1347989, C[12]);
    try std.testing.expectEqual(1538798, C[13]);
    try std.testing.expectEqual(1729607, C[14]);
    try std.testing.expectEqual(1569128, C[15]);
    try std.testing.expectEqual(1611881, C[16]);
    try std.testing.expectEqual(1862594, C[17]);
    try std.testing.expectEqual(2113307, C[18]);
    try std.testing.expectEqual(2364020, C[19]);

    const gamma = cf64.init(2, 2);
    const delta = cf64.init(3, 3);

    const D = try a.alloc(cf64, m * k);
    defer a.free(D);
    const E = try a.alloc(cf64, k * n);
    defer a.free(E);
    const F = try a.alloc(cf64, m * n);
    defer a.free(F);

    @memcpy(D.ptr, &[_]cf64{
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
    });
    @memcpy(E.ptr, &[_]cf64{
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
    @memcpy(F.ptr, &[_]cf64{
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
    });

    gemm(.row_major, .no_trans, .no_trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-184, F[0].re);
    try std.testing.expectEqual(190, F[0].im);
    try std.testing.expectEqual(-208, F[1].re);
    try std.testing.expectEqual(220, F[1].im);
    try std.testing.expectEqual(-232, F[2].re);
    try std.testing.expectEqual(250, F[2].im);
    try std.testing.expectEqual(-256, F[3].re);
    try std.testing.expectEqual(280, F[3].im);
    try std.testing.expectEqual(-280, F[4].re);
    try std.testing.expectEqual(310, F[4].im);
    try std.testing.expectEqual(-400, F[5].re);
    try std.testing.expectEqual(436, F[5].im);
    try std.testing.expectEqual(-460, F[6].re);
    try std.testing.expectEqual(502, F[6].im);
    try std.testing.expectEqual(-520, F[7].re);
    try std.testing.expectEqual(568, F[7].im);
    try std.testing.expectEqual(-580, F[8].re);
    try std.testing.expectEqual(634, F[8].im);
    try std.testing.expectEqual(-640, F[9].re);
    try std.testing.expectEqual(700, F[9].im);
    try std.testing.expectEqual(-616, F[10].re);
    try std.testing.expectEqual(682, F[10].im);
    try std.testing.expectEqual(-712, F[11].re);
    try std.testing.expectEqual(784, F[11].im);
    try std.testing.expectEqual(-808, F[12].re);
    try std.testing.expectEqual(886, F[12].im);
    try std.testing.expectEqual(-904, F[13].re);
    try std.testing.expectEqual(988, F[13].im);
    try std.testing.expectEqual(-1000, F[14].re);
    try std.testing.expectEqual(1090, F[14].im);
    try std.testing.expectEqual(-832, F[15].re);
    try std.testing.expectEqual(928, F[15].im);
    try std.testing.expectEqual(-964, F[16].re);
    try std.testing.expectEqual(1066, F[16].im);
    try std.testing.expectEqual(-1096, F[17].re);
    try std.testing.expectEqual(1204, F[17].im);
    try std.testing.expectEqual(-1228, F[18].re);
    try std.testing.expectEqual(1342, F[18].im);
    try std.testing.expectEqual(-1360, F[19].re);
    try std.testing.expectEqual(1480, F[19].im);

    gemm(.col_major, .no_trans, .no_trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-1274, F[0].re);
    try std.testing.expectEqual(170, F[0].im);
    try std.testing.expectEqual(-1460, F[1].re);
    try std.testing.expectEqual(212, F[1].im);
    try std.testing.expectEqual(-1646, F[2].re);
    try std.testing.expectEqual(254, F[2].im);
    try std.testing.expectEqual(-1832, F[3].re);
    try std.testing.expectEqual(296, F[3].im);
    try std.testing.expectEqual(-2102, F[4].re);
    try std.testing.expectEqual(422, F[4].im);
    try std.testing.expectEqual(-2900, F[5].re);
    try std.testing.expectEqual(500, F[5].im);
    try std.testing.expectEqual(-3338, F[6].re);
    try std.testing.expectEqual(578, F[6].im);
    try std.testing.expectEqual(-3776, F[7].re);
    try std.testing.expectEqual(656, F[7].im);
    try std.testing.expectEqual(-4154, F[8].re);
    try std.testing.expectEqual(674, F[8].im);
    try std.testing.expectEqual(-4628, F[9].re);
    try std.testing.expectEqual(788, F[9].im);
    try std.testing.expectEqual(-4598, F[10].re);
    try std.testing.expectEqual(902, F[10].im);
    try std.testing.expectEqual(-5288, F[11].re);
    try std.testing.expectEqual(1016, F[11].im);
    try std.testing.expectEqual(-5774, F[12].re);
    try std.testing.expectEqual(926, F[12].im);
    try std.testing.expectEqual(-6500, F[13].re);
    try std.testing.expectEqual(1076, F[13].im);
    try std.testing.expectEqual(-7226, F[14].re);
    try std.testing.expectEqual(1226, F[14].im);
    try std.testing.expectEqual(-6368, F[15].re);
    try std.testing.expectEqual(1376, F[15].im);
    try std.testing.expectEqual(-6962, F[16].re);
    try std.testing.expectEqual(1178, F[16].im);
    try std.testing.expectEqual(-7940, F[17].re);
    try std.testing.expectEqual(1364, F[17].im);
    try std.testing.expectEqual(-8918, F[18].re);
    try std.testing.expectEqual(1550, F[18].im);
    try std.testing.expectEqual(-9896, F[19].re);
    try std.testing.expectEqual(1736, F[19].im);

    gemm(.row_major, .no_trans, .conj_no_trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-4148, F[0].re);
    try std.testing.expectEqual(-3128, F[0].im);
    try std.testing.expectEqual(-4808, F[1].re);
    try std.testing.expectEqual(-3536, F[1].im);
    try std.testing.expectEqual(-5468, F[2].re);
    try std.testing.expectEqual(-3944, F[2].im);
    try std.testing.expectEqual(-6128, F[3].re);
    try std.testing.expectEqual(-4352, F[3].im);
    try std.testing.expectEqual(-7292, F[4].re);
    try std.testing.expectEqual(-4760, F[4].im);
    try std.testing.expectEqual(-9800, F[5].re);
    try std.testing.expectEqual(-6800, F[5].im);
    try std.testing.expectEqual(-11288, F[6].re);
    try std.testing.expectEqual(-7820, F[6].im);
    try std.testing.expectEqual(-12776, F[7].re);
    try std.testing.expectEqual(-8840, F[7].im);
    try std.testing.expectEqual(-13904, F[8].re);
    try std.testing.expectEqual(-9860, F[8].im);
    try std.testing.expectEqual(-15608, F[9].re);
    try std.testing.expectEqual(-10880, F[9].im);
    try std.testing.expectEqual(-15884, F[10].re);
    try std.testing.expectEqual(-10472, F[10].im);
    try std.testing.expectEqual(-18200, F[11].re);
    try std.testing.expectEqual(-12104, F[11].im);
    try std.testing.expectEqual(-19292, F[12].re);
    try std.testing.expectEqual(-13736, F[12].im);
    try std.testing.expectEqual(-21824, F[13].re);
    try std.testing.expectEqual(-15368, F[13].im);
    try std.testing.expectEqual(-24356, F[14].re);
    try std.testing.expectEqual(-17000, F[14].im);
    try std.testing.expectEqual(-22400, F[15].re);
    try std.testing.expectEqual(-14144, F[15].im);
    try std.testing.expectEqual(-23456, F[16].re);
    try std.testing.expectEqual(-16388, F[16].im);
    try std.testing.expectEqual(-26816, F[17].re);
    try std.testing.expectEqual(-18632, F[17].im);
    try std.testing.expectEqual(-30176, F[18].re);
    try std.testing.expectEqual(-20876, F[18].im);
    try std.testing.expectEqual(-33536, F[19].re);
    try std.testing.expectEqual(-23120, F[19].im);

    gemm(.col_major, .no_trans, .conj_no_trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-2908, F[0].re);
    try std.testing.expectEqual(-21676, F[0].im);
    try std.testing.expectEqual(-3640, F[1].re);
    try std.testing.expectEqual(-24856, F[1].im);
    try std.testing.expectEqual(-4372, F[2].re);
    try std.testing.expectEqual(-28036, F[2].im);
    try std.testing.expectEqual(-5104, F[3].re);
    try std.testing.expectEqual(-31216, F[3].im);
    try std.testing.expectEqual(-7264, F[4].re);
    try std.testing.expectEqual(-35824, F[4].im);
    try std.testing.expectEqual(-8608, F[5].re);
    try std.testing.expectEqual(-49408, F[5].im);
    try std.testing.expectEqual(-9952, F[6].re);
    try std.testing.expectEqual(-56872, F[6].im);
    try std.testing.expectEqual(-11296, F[7].re);
    try std.testing.expectEqual(-64336, F[7].im);
    try std.testing.expectEqual(-11620, F[8].re);
    try std.testing.expectEqual(-70780, F[8].im);
    try std.testing.expectEqual(-13576, F[9].re);
    try std.testing.expectEqual(-78856, F[9].im);
    try std.testing.expectEqual(-15532, F[10].re);
    try std.testing.expectEqual(-78364, F[10].im);
    try std.testing.expectEqual(-17488, F[11].re);
    try std.testing.expectEqual(-90112, F[11].im);
    try std.testing.expectEqual(-15976, F[12].re);
    try std.testing.expectEqual(-98392, F[12].im);
    try std.testing.expectEqual(-18544, F[13].re);
    try std.testing.expectEqual(-110752, F[13].im);
    try std.testing.expectEqual(-21112, F[14].re);
    try std.testing.expectEqual(-123112, F[14].im);
    try std.testing.expectEqual(-23680, F[15].re);
    try std.testing.expectEqual(-108544, F[15].im);
    try std.testing.expectEqual(-20332, F[16].re);
    try std.testing.expectEqual(-118660, F[16].im);
    try std.testing.expectEqual(-23512, F[17].re);
    try std.testing.expectEqual(-135304, F[17].im);
    try std.testing.expectEqual(-26692, F[18].re);
    try std.testing.expectEqual(-151948, F[18].im);
    try std.testing.expectEqual(-29872, F[19].re);
    try std.testing.expectEqual(-168592, F[19].im);

    gemm(.row_major, .no_trans, .trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(56248, F[0].re);
    try std.testing.expectEqual(-73696, F[0].im);
    try std.testing.expectEqual(63520, F[1].re);
    try std.testing.expectEqual(-85360, F[1].im);
    try std.testing.expectEqual(70792, F[2].re);
    try std.testing.expectEqual(-97024, F[2].im);
    try std.testing.expectEqual(78064, F[3].re);
    try std.testing.expectEqual(-108688, F[3].im);
    try std.testing.expectEqual(85336, F[4].re);
    try std.testing.expectEqual(-128920, F[4].im);
    try std.testing.expectEqual(122272, F[5].re);
    try std.testing.expectEqual(-173920, F[5].im);
    try std.testing.expectEqual(140452, F[6].re);
    try std.testing.expectEqual(-200164, F[6].im);
    try std.testing.expectEqual(158632, F[7].re);
    try std.testing.expectEqual(-226408, F[7].im);
    try std.testing.expectEqual(176812, F[8].re);
    try std.testing.expectEqual(-246532, F[8].im);
    try std.testing.expectEqual(194992, F[9].re);
    try std.testing.expectEqual(-276448, F[9].im);
    try std.testing.expectEqual(188296, F[10].re);
    try std.testing.expectEqual(-281488, F[10].im);
    try std.testing.expectEqual(217384, F[11].re);
    try std.testing.expectEqual(-322312, F[11].im);
    try std.testing.expectEqual(246472, F[12].re);
    try std.testing.expectEqual(-342328, F[12].im);
    try std.testing.expectEqual(275560, F[13].re);
    try std.testing.expectEqual(-386824, F[13].im);
    try std.testing.expectEqual(304648, F[14].re);
    try std.testing.expectEqual(-431320, F[14].im);
    try std.testing.expectEqual(254320, F[15].re);
    try std.testing.expectEqual(-396400, F[15].im);
    try std.testing.expectEqual(294316, F[16].re);
    try std.testing.expectEqual(-416308, F[16].im);
    try std.testing.expectEqual(334312, F[17].re);
    try std.testing.expectEqual(-475384, F[17].im);
    try std.testing.expectEqual(374308, F[18].re);
    try std.testing.expectEqual(-534460, F[18].im);
    try std.testing.expectEqual(414304, F[19].re);
    try std.testing.expectEqual(-593536, F[19].im);

    gemm(.col_major, .no_trans, .trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(389312, F[0].re);
    try std.testing.expectEqual(-51824, F[0].im);
    try std.testing.expectEqual(446048, F[1].re);
    try std.testing.expectEqual(-64928, F[1].im);
    try std.testing.expectEqual(502784, F[2].re);
    try std.testing.expectEqual(-78032, F[2].im);
    try std.testing.expectEqual(559520, F[3].re);
    try std.testing.expectEqual(-91136, F[3].im);
    try std.testing.expectEqual(642188, F[4].re);
    try std.testing.expectEqual(-130172, F[4].im);
    try std.testing.expectEqual(887912, F[5].re);
    try std.testing.expectEqual(-154280, F[5].im);
    try std.testing.expectEqual(1021100, F[6].re);
    try std.testing.expectEqual(-178388, F[6].im);
    try std.testing.expectEqual(1154288, F[7].re);
    try std.testing.expectEqual(-202496, F[7].im);
    try std.testing.expectEqual(1269392, F[8].re);
    try std.testing.expectEqual(-208520, F[8].im);
    try std.testing.expectEqual(1413584, F[9].re);
    try std.testing.expectEqual(-243632, F[9].im);
    try std.testing.expectEqual(1408520, F[10].re);
    try std.testing.expectEqual(-278744, F[10].im);
    try std.testing.expectEqual(1618160, F[11].re);
    try std.testing.expectEqual(-313856, F[11].im);
    try std.testing.expectEqual(1765700, F[12].re);
    try std.testing.expectEqual(-286868, F[12].im);
    try std.testing.expectEqual(1986344, F[13].re);
    try std.testing.expectEqual(-332984, F[13].im);
    try std.testing.expectEqual(2206988, F[14].re);
    try std.testing.expectEqual(-379100, F[14].im);
    try std.testing.expectEqual(1951136, F[15].re);
    try std.testing.expectEqual(-425216, F[15].im);
    try std.testing.expectEqual(2131112, F[16].re);
    try std.testing.expectEqual(-365216, F[16].im);
    try std.testing.expectEqual(2428208, F[17].re);
    try std.testing.expectEqual(-422336, F[17].im);
    try std.testing.expectEqual(2725304, F[18].re);
    try std.testing.expectEqual(-479456, F[18].im);
    try std.testing.expectEqual(3022400, F[19].re);
    try std.testing.expectEqual(-536576, F[19].im);

    gemm(.row_major, .no_trans, .conj_trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(1323464, F[0].re);
    try std.testing.expectEqual(1012520, F[0].im);
    try std.testing.expectEqual(1533056, F[1].re);
    try std.testing.expectEqual(1143488, F[1].im);
    try std.testing.expectEqual(1742648, F[2].re);
    try std.testing.expectEqual(1274456, F[2].im);
    try std.testing.expectEqual(1952240, F[3].re);
    try std.testing.expectEqual(1405424, F[3].im);
    try std.testing.expectEqual(2317424, F[4].re);
    try std.testing.expectEqual(1536392, F[4].im);
    try std.testing.expectEqual(3126704, F[5].re);
    try std.testing.expectEqual(2201024, F[5].im);
    try std.testing.expectEqual(3598772, F[6].re);
    try std.testing.expectEqual(2528444, F[6].im);
    try std.testing.expectEqual(4070840, F[7].re);
    try std.testing.expectEqual(2855864, F[7].im);
    try std.testing.expectEqual(4434404, F[8].re);
    try std.testing.expectEqual(3183284, F[8].im);
    try std.testing.expectEqual(4972496, F[9].re);
    try std.testing.expectEqual(3510704, F[9].im);
    try std.testing.expectEqual(5061992, F[10].re);
    try std.testing.expectEqual(3389528, F[10].im);
    try std.testing.expectEqual(5796536, F[11].re);
    try std.testing.expectEqual(3913400, F[11].im);
    try std.testing.expectEqual(6158480, F[12].re);
    try std.testing.expectEqual(4437272, F[12].im);
    try std.testing.expectEqual(6959048, F[13].re);
    try std.testing.expectEqual(4961144, F[13].im);
    try std.testing.expectEqual(7759616, F[14].re);
    try std.testing.expectEqual(5485016, F[14].im);
    try std.testing.expectEqual(7129328, F[15].re);
    try std.testing.expectEqual(4578032, F[15].im);
    try std.testing.expectEqual(7489652, F[16].re);
    try std.testing.expectEqual(5298356, F[16].im);
    try std.testing.expectEqual(8552696, F[17].re);
    try std.testing.expectEqual(6018680, F[17].im);
    try std.testing.expectEqual(9615740, F[18].re);
    try std.testing.expectEqual(6739004, F[18].im);
    try std.testing.expectEqual(10678784, F[19].re);
    try std.testing.expectEqual(7459328, F[19].im);

    gemm(.col_major, .no_trans, .conj_trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(933352, F[0].re);
    try std.testing.expectEqual(7008472, F[0].im);
    try std.testing.expectEqual(1169296, F[1].re);
    try std.testing.expectEqual(8030224, F[1].im);
    try std.testing.expectEqual(1405240, F[2].re);
    try std.testing.expectEqual(9051976, F[2].im);
    try std.testing.expectEqual(1641184, F[3].re);
    try std.testing.expectEqual(10073728, F[3].im);
    try std.testing.expectEqual(2343676, F[4].re);
    try std.testing.expectEqual(11562028, F[4].im);
    try std.testing.expectEqual(2777704, F[5].re);
    try std.testing.expectEqual(15983848, F[5].im);
    try std.testing.expectEqual(3211732, F[6].re);
    try std.testing.expectEqual(18382396, F[6].im);
    try std.testing.expectEqual(3645760, F[7].re);
    try std.testing.expectEqual(20780944, F[7].im);
    try std.testing.expectEqual(3754000, F[8].re);
    try std.testing.expectEqual(22853704, F[8].im);
    try std.testing.expectEqual(4386112, F[9].re);
    try std.testing.expectEqual(25450336, F[9].im);
    try std.testing.expectEqual(5018224, F[10].re);
    try std.testing.expectEqual(25355392, F[10].im);
    try std.testing.expectEqual(5650336, F[11].re);
    try std.testing.expectEqual(29130736, F[11].im);
    try std.testing.expectEqual(5164324, F[12].re);
    try std.testing.expectEqual(31787956, F[12].im);
    try std.testing.expectEqual(5994520, F[13].re);
    try std.testing.expectEqual(35761384, F[13].im);
    try std.testing.expectEqual(6824716, F[14].re);
    try std.testing.expectEqual(39734812, F[14].im);
    try std.testing.expectEqual(7654912, F[15].re);
    try std.testing.expectEqual(35123104, F[15].im);
    try std.testing.expectEqual(6574648, F[16].re);
    try std.testing.expectEqual(38364784, F[16].im);
    try std.testing.expectEqual(7602928, F[17].re);
    try std.testing.expectEqual(43715008, F[17].im);
    try std.testing.expectEqual(8631208, F[18].re);
    try std.testing.expectEqual(49065232, F[18].im);
    try std.testing.expectEqual(9659488, F[19].re);
    try std.testing.expectEqual(54415456, F[19].im);

    @memcpy(F.ptr, &[_]cf64{
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
    });

    gemm(.row_major, .conj_no_trans, .no_trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(184, F[0].re);
    try std.testing.expectEqual(190, F[0].im);
    try std.testing.expectEqual(208, F[1].re);
    try std.testing.expectEqual(220, F[1].im);
    try std.testing.expectEqual(232, F[2].re);
    try std.testing.expectEqual(250, F[2].im);
    try std.testing.expectEqual(256, F[3].re);
    try std.testing.expectEqual(280, F[3].im);
    try std.testing.expectEqual(280, F[4].re);
    try std.testing.expectEqual(310, F[4].im);
    try std.testing.expectEqual(400, F[5].re);
    try std.testing.expectEqual(436, F[5].im);
    try std.testing.expectEqual(460, F[6].re);
    try std.testing.expectEqual(502, F[6].im);
    try std.testing.expectEqual(520, F[7].re);
    try std.testing.expectEqual(568, F[7].im);
    try std.testing.expectEqual(580, F[8].re);
    try std.testing.expectEqual(634, F[8].im);
    try std.testing.expectEqual(640, F[9].re);
    try std.testing.expectEqual(700, F[9].im);
    try std.testing.expectEqual(616, F[10].re);
    try std.testing.expectEqual(682, F[10].im);
    try std.testing.expectEqual(712, F[11].re);
    try std.testing.expectEqual(784, F[11].im);
    try std.testing.expectEqual(808, F[12].re);
    try std.testing.expectEqual(886, F[12].im);
    try std.testing.expectEqual(904, F[13].re);
    try std.testing.expectEqual(988, F[13].im);
    try std.testing.expectEqual(1000, F[14].re);
    try std.testing.expectEqual(1090, F[14].im);
    try std.testing.expectEqual(832, F[15].re);
    try std.testing.expectEqual(928, F[15].im);
    try std.testing.expectEqual(964, F[16].re);
    try std.testing.expectEqual(1066, F[16].im);
    try std.testing.expectEqual(1096, F[17].re);
    try std.testing.expectEqual(1204, F[17].im);
    try std.testing.expectEqual(1228, F[18].re);
    try std.testing.expectEqual(1342, F[18].im);
    try std.testing.expectEqual(1360, F[19].re);
    try std.testing.expectEqual(1480, F[19].im);

    gemm(.col_major, .conj_no_trans, .no_trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(134, F[0].re);
    try std.testing.expectEqual(1274, F[0].im);
    try std.testing.expectEqual(140, F[1].re);
    try std.testing.expectEqual(1460, F[1].im);
    try std.testing.expectEqual(146, F[2].re);
    try std.testing.expectEqual(1646, F[2].im);
    try std.testing.expectEqual(152, F[3].re);
    try std.testing.expectEqual(1832, F[3].im);
    try std.testing.expectEqual(242, F[4].re);
    try std.testing.expectEqual(2102, F[4].im);
    try std.testing.expectEqual(284, F[5].re);
    try std.testing.expectEqual(2900, F[5].im);
    try std.testing.expectEqual(326, F[6].re);
    try std.testing.expectEqual(3338, F[6].im);
    try std.testing.expectEqual(368, F[7].re);
    try std.testing.expectEqual(3776, F[7].im);
    try std.testing.expectEqual(350, F[8].re);
    try std.testing.expectEqual(4154, F[8].im);
    try std.testing.expectEqual(428, F[9].re);
    try std.testing.expectEqual(4628, F[9].im);
    try std.testing.expectEqual(506, F[10].re);
    try std.testing.expectEqual(4598, F[10].im);
    try std.testing.expectEqual(584, F[11].re);
    try std.testing.expectEqual(5288, F[11].im);
    try std.testing.expectEqual(458, F[12].re);
    try std.testing.expectEqual(5774, F[12].im);
    try std.testing.expectEqual(572, F[13].re);
    try std.testing.expectEqual(6500, F[13].im);
    try std.testing.expectEqual(686, F[14].re);
    try std.testing.expectEqual(7226, F[14].im);
    try std.testing.expectEqual(800, F[15].re);
    try std.testing.expectEqual(6368, F[15].im);
    try std.testing.expectEqual(566, F[16].re);
    try std.testing.expectEqual(6962, F[16].im);
    try std.testing.expectEqual(716, F[17].re);
    try std.testing.expectEqual(7940, F[17].im);
    try std.testing.expectEqual(866, F[18].re);
    try std.testing.expectEqual(8918, F[18].im);
    try std.testing.expectEqual(1016, F[19].re);
    try std.testing.expectEqual(9896, F[19].im);

    gemm(.row_major, .conj_no_trans, .conj_no_trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-3236, F[0].re);
    try std.testing.expectEqual(4040, F[0].im);
    try std.testing.expectEqual(-3752, F[1].re);
    try std.testing.expectEqual(4592, F[1].im);
    try std.testing.expectEqual(-4268, F[2].re);
    try std.testing.expectEqual(5144, F[2].im);
    try std.testing.expectEqual(-4784, F[3].re);
    try std.testing.expectEqual(5696, F[3].im);
    try std.testing.expectEqual(-5300, F[4].re);
    try std.testing.expectEqual(6752, F[4].im);
    try std.testing.expectEqual(-7448, F[5].re);
    try std.testing.expectEqual(9152, F[5].im);
    try std.testing.expectEqual(-8576, F[6].re);
    try std.testing.expectEqual(10532, F[6].im);
    try std.testing.expectEqual(-9704, F[7].re);
    try std.testing.expectEqual(11912, F[7].im);
    try std.testing.expectEqual(-10832, F[8].re);
    try std.testing.expectEqual(12932, F[8].im);
    try std.testing.expectEqual(-11960, F[9].re);
    try std.testing.expectEqual(14528, F[9].im);
    try std.testing.expectEqual(-11660, F[10].re);
    try std.testing.expectEqual(14696, F[10].im);
    try std.testing.expectEqual(-13400, F[11].re);
    try std.testing.expectEqual(16904, F[11].im);
    try std.testing.expectEqual(-15140, F[12].re);
    try std.testing.expectEqual(17888, F[12].im);
    try std.testing.expectEqual(-16880, F[13].re);
    try std.testing.expectEqual(20312, F[13].im);
    try std.testing.expectEqual(-18620, F[14].re);
    try std.testing.expectEqual(22736, F[14].im);
    try std.testing.expectEqual(-15872, F[15].re);
    try std.testing.expectEqual(20672, F[15].im);
    try std.testing.expectEqual(-18224, F[16].re);
    try std.testing.expectEqual(21620, F[16].im);
    try std.testing.expectEqual(-20576, F[17].re);
    try std.testing.expectEqual(24872, F[17].im);
    try std.testing.expectEqual(-22928, F[18].re);
    try std.testing.expectEqual(28124, F[18].im);
    try std.testing.expectEqual(-25280, F[19].re);
    try std.testing.expectEqual(31376, F[19].im);

    gemm(.col_major, .conj_no_trans, .conj_no_trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-21676, F[0].re);
    try std.testing.expectEqual(2260, F[0].im);
    try std.testing.expectEqual(-24856, F[1].re);
    try std.testing.expectEqual(2344, F[1].im);
    try std.testing.expectEqual(-28036, F[2].re);
    try std.testing.expectEqual(2428, F[2].im);
    try std.testing.expectEqual(-31216, F[3].re);
    try std.testing.expectEqual(2512, F[3].im);
    try std.testing.expectEqual(-35824, F[4].re);
    try std.testing.expectEqual(4024, F[4].im);
    try std.testing.expectEqual(-49408, F[5].re);
    try std.testing.expectEqual(4720, F[5].im);
    try std.testing.expectEqual(-56872, F[6].re);
    try std.testing.expectEqual(5416, F[6].im);
    try std.testing.expectEqual(-64336, F[7].re);
    try std.testing.expectEqual(6112, F[7].im);
    try std.testing.expectEqual(-70780, F[8].re);
    try std.testing.expectEqual(5788, F[8].im);
    try std.testing.expectEqual(-78856, F[9].re);
    try std.testing.expectEqual(7096, F[9].im);
    try std.testing.expectEqual(-78364, F[10].re);
    try std.testing.expectEqual(8404, F[10].im);
    try std.testing.expectEqual(-90112, F[11].re);
    try std.testing.expectEqual(9712, F[11].im);
    try std.testing.expectEqual(-98392, F[12].re);
    try std.testing.expectEqual(7552, F[12].im);
    try std.testing.expectEqual(-110752, F[13].re);
    try std.testing.expectEqual(9472, F[13].im);
    try std.testing.expectEqual(-123112, F[14].re);
    try std.testing.expectEqual(11392, F[14].im);
    try std.testing.expectEqual(-108544, F[15].re);
    try std.testing.expectEqual(13312, F[15].im);
    try std.testing.expectEqual(-118660, F[16].re);
    try std.testing.expectEqual(9316, F[16].im);
    try std.testing.expectEqual(-135304, F[17].re);
    try std.testing.expectEqual(11848, F[17].im);
    try std.testing.expectEqual(-151948, F[18].re);
    try std.testing.expectEqual(14380, F[18].im);
    try std.testing.expectEqual(-168592, F[19].re);
    try std.testing.expectEqual(16912, F[19].im);

    gemm(.row_major, .conj_no_trans, .trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-71752, F[0].re);
    try std.testing.expectEqual(-58192, F[0].im);
    try std.testing.expectEqual(-81472, F[1].re);
    try std.testing.expectEqual(-67408, F[1].im);
    try std.testing.expectEqual(-91192, F[2].re);
    try std.testing.expectEqual(-76624, F[2].im);
    try std.testing.expectEqual(-100912, F[3].re);
    try std.testing.expectEqual(-85840, F[3].im);
    try std.testing.expectEqual(-119200, F[4].re);
    try std.testing.expectEqual(-95056, F[4].im);
    try std.testing.expectEqual(-162256, F[5].re);
    try std.testing.expectEqual(-133936, F[5].im);
    try std.testing.expectEqual(-186556, F[6].re);
    try std.testing.expectEqual(-154060, F[6].im);
    try std.testing.expectEqual(-210856, F[7].re);
    try std.testing.expectEqual(-174184, F[7].im);
    try std.testing.expectEqual(-229036, F[8].re);
    try std.testing.expectEqual(-194308, F[8].im);
    try std.testing.expectEqual(-257008, F[9].re);
    try std.testing.expectEqual(-214432, F[9].im);
    try std.testing.expectEqual(-260104, F[10].re);
    try std.testing.expectEqual(-209680, F[10].im);
    try std.testing.expectEqual(-298984, F[11].re);
    try std.testing.expectEqual(-240712, F[11].im);
    try std.testing.expectEqual(-317056, F[12].re);
    try std.testing.expectEqual(-271744, F[12].im);
    try std.testing.expectEqual(-359608, F[13].re);
    try std.testing.expectEqual(-302776, F[13].im);
    try std.testing.expectEqual(-402160, F[14].re);
    try std.testing.expectEqual(-333808, F[14].im);
    try std.testing.expectEqual(-365296, F[15].re);
    try std.testing.expectEqual(-285424, F[15].im);
    try std.testing.expectEqual(-383260, F[16].re);
    try std.testing.expectEqual(-327364, F[16].im);
    try std.testing.expectEqual(-440392, F[17].re);
    try std.testing.expectEqual(-369304, F[17].im);
    try std.testing.expectEqual(-497524, F[18].re);
    try std.testing.expectEqual(-411244, F[18].im);
    try std.testing.expectEqual(-554656, F[19].re);
    try std.testing.expectEqual(-453184, F[19].im);

    gemm(.col_major, .conj_no_trans, .trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-40160, F[0].re);
    try std.testing.expectEqual(-389312, F[0].im);
    try std.testing.expectEqual(-41600, F[1].re);
    try std.testing.expectEqual(-446048, F[1].im);
    try std.testing.expectEqual(-43040, F[2].re);
    try std.testing.expectEqual(-502784, F[2].im);
    try std.testing.expectEqual(-44480, F[3].re);
    try std.testing.expectEqual(-559520, F[3].im);
    try std.testing.expectEqual(-71852, F[4].re);
    try std.testing.expectEqual(-642188, F[4].im);
    try std.testing.expectEqual(-84296, F[5].re);
    try std.testing.expectEqual(-887912, F[5].im);
    try std.testing.expectEqual(-96740, F[6].re);
    try std.testing.expectEqual(-1021100, F[6].im);
    try std.testing.expectEqual(-109184, F[7].re);
    try std.testing.expectEqual(-1154288, F[7].im);
    try std.testing.expectEqual(-103544, F[8].re);
    try std.testing.expectEqual(-1269392, F[8].im);
    try std.testing.expectEqual(-126992, F[9].re);
    try std.testing.expectEqual(-1413584, F[9].im);
    try std.testing.expectEqual(-150440, F[10].re);
    try std.testing.expectEqual(-1408520, F[10].im);
    try std.testing.expectEqual(-173888, F[11].re);
    try std.testing.expectEqual(-1618160, F[11].im);
    try std.testing.expectEqual(-135236, F[12].re);
    try std.testing.expectEqual(-1765700, F[12].im);
    try std.testing.expectEqual(-169688, F[13].re);
    try std.testing.expectEqual(-1986344, F[13].im);
    try std.testing.expectEqual(-204140, F[14].re);
    try std.testing.expectEqual(-2206988, F[14].im);
    try std.testing.expectEqual(-238592, F[15].re);
    try std.testing.expectEqual(-1951136, F[15].im);
    try std.testing.expectEqual(-166928, F[16].re);
    try std.testing.expectEqual(-2131112, F[16].im);
    try std.testing.expectEqual(-212384, F[17].re);
    try std.testing.expectEqual(-2428208, F[17].im);
    try std.testing.expectEqual(-257840, F[18].re);
    try std.testing.expectEqual(-2725304, F[18].im);
    try std.testing.expectEqual(-303296, F[19].re);
    try std.testing.expectEqual(-3022400, F[19].im);

    gemm(.row_major, .conj_no_trans, .conj_trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(1047512, F[0].re);
    try std.testing.expectEqual(-1288472, F[0].im);
    try std.testing.expectEqual(1213472, F[1].re);
    try std.testing.expectEqual(-1463072, F[1].im);
    try std.testing.expectEqual(1379432, F[2].re);
    try std.testing.expectEqual(-1637672, F[2].im);
    try std.testing.expectEqual(1545392, F[3].re);
    try std.testing.expectEqual(-1812272, F[3].im);
    try std.testing.expectEqual(1711352, F[4].re);
    try std.testing.expectEqual(-2142464, F[4].im);
    try std.testing.expectEqual(2410976, F[5].re);
    try std.testing.expectEqual(-2916752, F[5].im);
    try std.testing.expectEqual(2773388, F[6].re);
    try std.testing.expectEqual(-3353828, F[6].im);
    try std.testing.expectEqual(3135800, F[7].re);
    try std.testing.expectEqual(-3790904, F[7].im);
    try std.testing.expectEqual(3498212, F[8].re);
    try std.testing.expectEqual(-4119476, F[8].im);
    try std.testing.expectEqual(3860624, F[9].re);
    try std.testing.expectEqual(-4622576, F[9].im);
    try std.testing.expectEqual(3774440, F[10].re);
    try std.testing.expectEqual(-4677080, F[10].im);
    try std.testing.expectEqual(4333304, F[11].re);
    try std.testing.expectEqual(-5376632, F[11].im);
    try std.testing.expectEqual(4892168, F[12].re);
    try std.testing.expectEqual(-5703584, F[12].im);
    try std.testing.expectEqual(5451032, F[13].re);
    try std.testing.expectEqual(-6469160, F[13].im);
    try std.testing.expectEqual(6009896, F[14].re);
    try std.testing.expectEqual(-7234736, F[14].im);
    try std.testing.expectEqual(5137904, F[15].re);
    try std.testing.expectEqual(-6569456, F[15].im);
    try std.testing.expectEqual(5893220, F[16].re);
    try std.testing.expectEqual(-6894788, F[16].im);
    try std.testing.expectEqual(6648536, F[17].re);
    try std.testing.expectEqual(-7922840, F[17].im);
    try std.testing.expectEqual(7403852, F[18].re);
    try std.testing.expectEqual(-8950892, F[18].im);
    try std.testing.expectEqual(8159168, F[19].re);
    try std.testing.expectEqual(-9978944, F[19].im);

    gemm(.col_major, .conj_no_trans, .conj_trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(7008472, F[0].re);
    try std.testing.expectEqual(-723400, F[0].im);
    try std.testing.expectEqual(8030224, F[1].re);
    try std.testing.expectEqual(-749392, F[1].im);
    try std.testing.expectEqual(9051976, F[2].re);
    try std.testing.expectEqual(-775384, F[2].im);
    try std.testing.expectEqual(10073728, F[3].re);
    try std.testing.expectEqual(-801376, F[3].im);
    try std.testing.expectEqual(11562028, F[4].re);
    try std.testing.expectEqual(-1293916, F[4].im);
    try std.testing.expectEqual(15983848, F[5].re);
    try std.testing.expectEqual(-1517992, F[5].im);
    try std.testing.expectEqual(18382396, F[6].re);
    try std.testing.expectEqual(-1742068, F[6].im);
    try std.testing.expectEqual(20780944, F[7].re);
    try std.testing.expectEqual(-1966144, F[7].im);
    try std.testing.expectEqual(22853704, F[8].re);
    try std.testing.expectEqual(-1864432, F[8].im);
    try std.testing.expectEqual(25450336, F[9].re);
    try std.testing.expectEqual(-2286592, F[9].im);
    try std.testing.expectEqual(25355392, F[10].re);
    try std.testing.expectEqual(-2708752, F[10].im);
    try std.testing.expectEqual(29130736, F[11].re);
    try std.testing.expectEqual(-3130912, F[11].im);
    try std.testing.expectEqual(31787956, F[12].re);
    try std.testing.expectEqual(-2434948, F[12].im);
    try std.testing.expectEqual(35761384, F[13].re);
    try std.testing.expectEqual(-3055192, F[13].im);
    try std.testing.expectEqual(39734812, F[14].re);
    try std.testing.expectEqual(-3675436, F[14].im);
    try std.testing.expectEqual(35123104, F[15].re);
    try std.testing.expectEqual(-4295680, F[15].im);
    try std.testing.expectEqual(38364784, F[16].re);
    try std.testing.expectEqual(-3005464, F[16].im);
    try std.testing.expectEqual(43715008, F[17].re);
    try std.testing.expectEqual(-3823792, F[17].im);
    try std.testing.expectEqual(49065232, F[18].re);
    try std.testing.expectEqual(-4642120, F[18].im);
    try std.testing.expectEqual(54415456, F[19].re);
    try std.testing.expectEqual(-5460448, F[19].im);

    @memcpy(F.ptr, &[_]cf64{
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
    });

    gemm(.row_major, .trans, .no_trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-520, F[0].re);
    try std.testing.expectEqual(526, F[0].im);
    try std.testing.expectEqual(-580, F[1].re);
    try std.testing.expectEqual(592, F[1].im);
    try std.testing.expectEqual(-640, F[2].re);
    try std.testing.expectEqual(658, F[2].im);
    try std.testing.expectEqual(-700, F[3].re);
    try std.testing.expectEqual(724, F[3].im);
    try std.testing.expectEqual(-760, F[4].re);
    try std.testing.expectEqual(790, F[4].im);
    try std.testing.expectEqual(-592, F[5].re);
    try std.testing.expectEqual(628, F[5].im);
    try std.testing.expectEqual(-664, F[6].re);
    try std.testing.expectEqual(706, F[6].im);
    try std.testing.expectEqual(-736, F[7].re);
    try std.testing.expectEqual(784, F[7].im);
    try std.testing.expectEqual(-808, F[8].re);
    try std.testing.expectEqual(862, F[8].im);
    try std.testing.expectEqual(-880, F[9].re);
    try std.testing.expectEqual(940, F[9].im);
    try std.testing.expectEqual(-664, F[10].re);
    try std.testing.expectEqual(730, F[10].im);
    try std.testing.expectEqual(-748, F[11].re);
    try std.testing.expectEqual(820, F[11].im);
    try std.testing.expectEqual(-832, F[12].re);
    try std.testing.expectEqual(910, F[12].im);
    try std.testing.expectEqual(-916, F[13].re);
    try std.testing.expectEqual(1000, F[13].im);
    try std.testing.expectEqual(-1000, F[14].re);
    try std.testing.expectEqual(1090, F[14].im);
    try std.testing.expectEqual(-736, F[15].re);
    try std.testing.expectEqual(832, F[15].im);
    try std.testing.expectEqual(-832, F[16].re);
    try std.testing.expectEqual(934, F[16].im);
    try std.testing.expectEqual(-928, F[17].re);
    try std.testing.expectEqual(1036, F[17].im);
    try std.testing.expectEqual(-1024, F[18].re);
    try std.testing.expectEqual(1138, F[18].im);
    try std.testing.expectEqual(-1120, F[19].re);
    try std.testing.expectEqual(1240, F[19].im);

    gemm(.col_major, .trans, .no_trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-3194, F[0].re);
    try std.testing.expectEqual(74, F[0].im);
    try std.testing.expectEqual(-3644, F[1].re);
    try std.testing.expectEqual(164, F[1].im);
    try std.testing.expectEqual(-4094, F[2].re);
    try std.testing.expectEqual(254, F[2].im);
    try std.testing.expectEqual(-4544, F[3].re);
    try std.testing.expectEqual(344, F[3].im);
    try std.testing.expectEqual(-4778, F[4].re);
    try std.testing.expectEqual(218, F[4].im);
    try std.testing.expectEqual(-3968, F[5].re);
    try std.testing.expectEqual(416, F[5].im);
    try std.testing.expectEqual(-4598, F[6].re);
    try std.testing.expectEqual(614, F[6].im);
    try std.testing.expectEqual(-5228, F[7].re);
    try std.testing.expectEqual(812, F[7].im);
    try std.testing.expectEqual(-5210, F[8].re);
    try std.testing.expectEqual(362, F[8].im);
    try std.testing.expectEqual(-5948, F[9].re);
    try std.testing.expectEqual(668, F[9].im);
    try std.testing.expectEqual(-4958, F[10].re);
    try std.testing.expectEqual(974, F[10].im);
    try std.testing.expectEqual(-5768, F[11].re);
    try std.testing.expectEqual(1280, F[11].im);
    try std.testing.expectEqual(-5498, F[12].re);
    try std.testing.expectEqual(506, F[12].im);
    try std.testing.expectEqual(-6416, F[13].re);
    try std.testing.expectEqual(920, F[13].im);
    try std.testing.expectEqual(-7334, F[14].re);
    try std.testing.expectEqual(1334, F[14].im);
    try std.testing.expectEqual(-6164, F[15].re);
    try std.testing.expectEqual(1748, F[15].im);
    try std.testing.expectEqual(-5642, F[16].re);
    try std.testing.expectEqual(650, F[16].im);
    try std.testing.expectEqual(-6740, F[17].re);
    try std.testing.expectEqual(1172, F[17].im);
    try std.testing.expectEqual(-7838, F[18].re);
    try std.testing.expectEqual(1694, F[18].im);
    try std.testing.expectEqual(-8936, F[19].re);
    try std.testing.expectEqual(2216, F[19].im);

    gemm(.row_major, .trans, .conj_no_trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-9284, F[0].re);
    try std.testing.expectEqual(-8840, F[0].im);
    try std.testing.expectEqual(-10844, F[1].re);
    try std.testing.expectEqual(-9860, F[1].im);
    try std.testing.expectEqual(-12404, F[2].re);
    try std.testing.expectEqual(-10880, F[2].im);
    try std.testing.expectEqual(-13964, F[3].re);
    try std.testing.expectEqual(-11900, F[3].im);
    try std.testing.expectEqual(-14228, F[4].re);
    try std.testing.expectEqual(-12920, F[4].im);
    try std.testing.expectEqual(-12560, F[5].re);
    try std.testing.expectEqual(-10064, F[5].im);
    try std.testing.expectEqual(-14972, F[6].re);
    try std.testing.expectEqual(-11288, F[6].im);
    try std.testing.expectEqual(-17384, F[7].re);
    try std.testing.expectEqual(-12512, F[7].im);
    try std.testing.expectEqual(-15908, F[8].re);
    try std.testing.expectEqual(-13736, F[8].im);
    try std.testing.expectEqual(-18968, F[9].re);
    try std.testing.expectEqual(-14960, F[9].im);
    try std.testing.expectEqual(-17132, F[10].re);
    try std.testing.expectEqual(-11288, F[10].im);
    try std.testing.expectEqual(-20396, F[11].re);
    try std.testing.expectEqual(-12716, F[11].im);
    try std.testing.expectEqual(-17180, F[12].re);
    try std.testing.expectEqual(-14144, F[12].im);
    try std.testing.expectEqual(-21092, F[13].re);
    try std.testing.expectEqual(-15572, F[13].im);
    try std.testing.expectEqual(-25004, F[14].re);
    try std.testing.expectEqual(-17000, F[14].im);
    try std.testing.expectEqual(-23000, F[15].re);
    try std.testing.expectEqual(-12512, F[15].im);
    try std.testing.expectEqual(-18044, F[16].re);
    try std.testing.expectEqual(-14144, F[16].im);
    try std.testing.expectEqual(-22808, F[17].re);
    try std.testing.expectEqual(-15776, F[17].im);
    try std.testing.expectEqual(-27572, F[18].re);
    try std.testing.expectEqual(-17408, F[18].im);
    try std.testing.expectEqual(-32336, F[19].re);
    try std.testing.expectEqual(-19040, F[19].im);

    gemm(.col_major, .trans, .conj_no_trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-1276, F[0].re);
    try std.testing.expectEqual(-54316, F[0].im);
    try std.testing.expectEqual(-2824, F[1].re);
    try std.testing.expectEqual(-61984, F[1].im);
    try std.testing.expectEqual(-4372, F[2].re);
    try std.testing.expectEqual(-69652, F[2].im);
    try std.testing.expectEqual(-5920, F[3].re);
    try std.testing.expectEqual(-77320, F[3].im);
    try std.testing.expectEqual(-3796, F[4].re);
    try std.testing.expectEqual(-81316, F[4].im);
    try std.testing.expectEqual(-7180, F[5].re);
    try std.testing.expectEqual(-67564, F[5].im);
    try std.testing.expectEqual(-10564, F[6].re);
    try std.testing.expectEqual(-78292, F[6].im);
    try std.testing.expectEqual(-13948, F[7].re);
    try std.testing.expectEqual(-89020, F[7].im);
    try std.testing.expectEqual(-6316, F[8].re);
    try std.testing.expectEqual(-88732, F[8].im);
    try std.testing.expectEqual(-11536, F[9].re);
    try std.testing.expectEqual(-101296, F[9].im);
    try std.testing.expectEqual(-16756, F[10].re);
    try std.testing.expectEqual(-84484, F[10].im);
    try std.testing.expectEqual(-21976, F[11].re);
    try std.testing.expectEqual(-98272, F[11].im);
    try std.testing.expectEqual(-8836, F[12].re);
    try std.testing.expectEqual(-93700, F[12].im);
    try std.testing.expectEqual(-15892, F[13].re);
    try std.testing.expectEqual(-109324, F[13].im);
    try std.testing.expectEqual(-22948, F[14].re);
    try std.testing.expectEqual(-124948, F[14].im);
    try std.testing.expectEqual(-30004, F[15].re);
    try std.testing.expectEqual(-105076, F[15].im);
    try std.testing.expectEqual(-11356, F[16].re);
    try std.testing.expectEqual(-96220, F[16].im);
    try std.testing.expectEqual(-20248, F[17].re);
    try std.testing.expectEqual(-114904, F[17].im);
    try std.testing.expectEqual(-29140, F[18].re);
    try std.testing.expectEqual(-133588, F[18].im);
    try std.testing.expectEqual(-38032, F[19].re);
    try std.testing.expectEqual(-152272, F[19].im);

    gemm(.row_major, .trans, .trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(158968, F[0].re);
    try std.testing.expectEqual(-166624, F[0].im);
    try std.testing.expectEqual(177148, F[1].re);
    try std.testing.expectEqual(-194092, F[1].im);
    try std.testing.expectEqual(195328, F[2].re);
    try std.testing.expectEqual(-221560, F[2].im);
    try std.testing.expectEqual(213508, F[3].re);
    try std.testing.expectEqual(-249028, F[3].im);
    try std.testing.expectEqual(231688, F[4].re);
    try std.testing.expectEqual(-254464, F[4].im);
    try std.testing.expectEqual(180976, F[5].re);
    try std.testing.expectEqual(-224056, F[5].im);
    try std.testing.expectEqual(202792, F[6].re);
    try std.testing.expectEqual(-266176, F[6].im);
    try std.testing.expectEqual(224608, F[7].re);
    try std.testing.expectEqual(-308296, F[7].im);
    try std.testing.expectEqual(246424, F[8].re);
    try std.testing.expectEqual(-284320, F[8].im);
    try std.testing.expectEqual(268240, F[9].re);
    try std.testing.expectEqual(-337456, F[9].im);
    try std.testing.expectEqual(202984, F[10].re);
    try std.testing.expectEqual(-303520, F[10].im);
    try std.testing.expectEqual(228436, F[11].re);
    try std.testing.expectEqual(-360292, F[11].im);
    try std.testing.expectEqual(253888, F[12].re);
    try std.testing.expectEqual(-306904, F[12].im);
    try std.testing.expectEqual(279340, F[13].re);
    try std.testing.expectEqual(-374692, F[13].im);
    try std.testing.expectEqual(304792, F[14].re);
    try std.testing.expectEqual(-442480, F[14].im);
    try std.testing.expectEqual(224992, F[15].re);
    try std.testing.expectEqual(-405016, F[15].im);
    try std.testing.expectEqual(254080, F[16].re);
    try std.testing.expectEqual(-322216, F[16].im);
    try std.testing.expectEqual(283168, F[17].re);
    try std.testing.expectEqual(-404656, F[17].im);
    try std.testing.expectEqual(312256, F[18].re);
    try std.testing.expectEqual(-487096, F[18].im);
    try std.testing.expectEqual(341344, F[19].re);
    try std.testing.expectEqual(-569536, F[19].im);

    gemm(.col_major, .trans, .trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(976592, F[0].re);
    try std.testing.expectEqual(-22784, F[0].im);
    try std.testing.expectEqual(1113320, F[1].re);
    try std.testing.expectEqual(-50432, F[1].im);
    try std.testing.expectEqual(1250048, F[2].re);
    try std.testing.expectEqual(-78080, F[2].im);
    try std.testing.expectEqual(1386776, F[3].re);
    try std.testing.expectEqual(-105728, F[3].im);
    try std.testing.expectEqual(1458248, F[4].re);
    try std.testing.expectEqual(-68120, F[4].im);
    try std.testing.expectEqual(1214636, F[5].re);
    try std.testing.expectEqual(-128780, F[5].im);
    try std.testing.expectEqual(1406192, F[6].re);
    try std.testing.expectEqual(-189440, F[6].im);
    try std.testing.expectEqual(1597748, F[7].re);
    try std.testing.expectEqual(-250100, F[7].im);
    try std.testing.expectEqual(1592000, F[8].re);
    try std.testing.expectEqual(-113456, F[8].im);
    try std.testing.expectEqual(1816568, F[9].re);
    try std.testing.expectEqual(-207128, F[9].im);
    try std.testing.expectEqual(1518704, F[10].re);
    try std.testing.expectEqual(-300800, F[10].im);
    try std.testing.expectEqual(1765088, F[11].re);
    try std.testing.expectEqual(-394472, F[11].im);
    try std.testing.expectEqual(1682120, F[12].re);
    try std.testing.expectEqual(-158792, F[12].im);
    try std.testing.expectEqual(1961516, F[13].re);
    try std.testing.expectEqual(-285476, F[13].im);
    try std.testing.expectEqual(2240912, F[14].re);
    try std.testing.expectEqual(-412160, F[14].im);
    try std.testing.expectEqual(1888796, F[15].re);
    try std.testing.expectEqual(-538844, F[15].im);
    try std.testing.expectEqual(1728608, F[16].re);
    try std.testing.expectEqual(-204128, F[16].im);
    try std.testing.expectEqual(2062832, F[17].re);
    try std.testing.expectEqual(-363824, F[17].im);
    try std.testing.expectEqual(2397056, F[18].re);
    try std.testing.expectEqual(-523520, F[18].im);
    try std.testing.expectEqual(2731280, F[19].re);
    try std.testing.expectEqual(-683216, F[19].im);

    gemm(.row_major, .trans, .conj_trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(2998280, F[0].re);
    try std.testing.expectEqual(2861576, F[0].im);
    try std.testing.expectEqual(3491588, F[1].re);
    try std.testing.expectEqual(3188996, F[1].im);
    try std.testing.expectEqual(3984896, F[2].re);
    try std.testing.expectEqual(3516416, F[2].im);
    try std.testing.expectEqual(4478204, F[3].re);
    try std.testing.expectEqual(3843836, F[3].im);
    try std.testing.expectEqual(4579976, F[4].re);
    try std.testing.expectEqual(4171256, F[4].im);
    try std.testing.expectEqual(4030424, F[5].re);
    try std.testing.expectEqual(3257744, F[5].im);
    try std.testing.expectEqual(4787288, F[6].re);
    try std.testing.expectEqual(3650648, F[6].im);
    try std.testing.expectEqual(5544152, F[7].re);
    try std.testing.expectEqual(4043552, F[7].im);
    try std.testing.expectEqual(5117192, F[8].re);
    try std.testing.expectEqual(4436456, F[8].im);
    try std.testing.expectEqual(6072128, F[9].re);
    try std.testing.expectEqual(4829360, F[9].im);
    try std.testing.expectEqual(5458712, F[10].re);
    try std.testing.expectEqual(3653912, F[10].im);
    try std.testing.expectEqual(6479132, F[11].re);
    try std.testing.expectEqual(4112300, F[11].im);
    try std.testing.expectEqual(5523440, F[12].re);
    try std.testing.expectEqual(4570688, F[12].im);
    try std.testing.expectEqual(6741932, F[13].re);
    try std.testing.expectEqual(5029076, F[13].im);
    try std.testing.expectEqual(7960424, F[14].re);
    try std.testing.expectEqual(5487464, F[14].im);
    try std.testing.expectEqual(7283144, F[15].re);
    try std.testing.expectEqual(4050080, F[15].im);
    try std.testing.expectEqual(5798720, F[16].re);
    try std.testing.expectEqual(4573952, F[16].im);
    try std.testing.expectEqual(7280768, F[17].re);
    try std.testing.expectEqual(5097824, F[17].im);
    try std.testing.expectEqual(8762816, F[18].re);
    try std.testing.expectEqual(5621696, F[18].im);
    try std.testing.expectEqual(10244864, F[19].re);
    try std.testing.expectEqual(6145568, F[19].im);

    gemm(.col_major, .trans, .conj_trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(410296, F[0].re);
    try std.testing.expectEqual(17579752, F[0].im);
    try std.testing.expectEqual(908176, F[1].re);
    try std.testing.expectEqual(20042152, F[1].im);
    try std.testing.expectEqual(1406056, F[2].re);
    try std.testing.expectEqual(22504552, F[2].im);
    try std.testing.expectEqual(1903936, F[3].re);
    try std.testing.expectEqual(24966952, F[3].im);
    try std.testing.expectEqual(1226368, F[4].re);
    try std.testing.expectEqual(26253904, F[4].im);
    try std.testing.expectEqual(2318500, F[5].re);
    try std.testing.expectEqual(21864964, F[5].im);
    try std.testing.expectEqual(3410632, F[6].re);
    try std.testing.expectEqual(25314520, F[6].im);
    try std.testing.expectEqual(4502764, F[7].re);
    try std.testing.expectEqual(28764076, F[7].im);
    try std.testing.expectEqual(2042440, F[8].re);
    try std.testing.expectEqual(28661176, F[8].im);
    try std.testing.expectEqual(3728824, F[9].re);
    try std.testing.expectEqual(32704984, F[9].im);
    try std.testing.expectEqual(5415208, F[10].re);
    try std.testing.expectEqual(27338680, F[10].im);
    try std.testing.expectEqual(7101592, F[11].re);
    try std.testing.expectEqual(31775392, F[11].im);
    try std.testing.expectEqual(2858512, F[12].re);
    try std.testing.expectEqual(30282640, F[12].im);
    try std.testing.expectEqual(5139148, F[13].re);
    try std.testing.expectEqual(35313604, F[13].im);
    try std.testing.expectEqual(7419784, F[14].re);
    try std.testing.expectEqual(40344568, F[14].im);
    try std.testing.expectEqual(9700420, F[15].re);
    try std.testing.expectEqual(34000900, F[15].im);
    try std.testing.expectEqual(3674584, F[16].re);
    try std.testing.expectEqual(31118296, F[16].im);
    try std.testing.expectEqual(6549472, F[17].re);
    try std.testing.expectEqual(37136416, F[17].im);
    try std.testing.expectEqual(9424360, F[18].re);
    try std.testing.expectEqual(43154536, F[18].im);
    try std.testing.expectEqual(12299248, F[19].re);
    try std.testing.expectEqual(49172656, F[19].im);

    @memcpy(F.ptr, &[_]cf64{
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
    });

    gemm(.row_major, .conj_trans, .no_trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(520, F[0].re);
    try std.testing.expectEqual(526, F[0].im);
    try std.testing.expectEqual(580, F[1].re);
    try std.testing.expectEqual(592, F[1].im);
    try std.testing.expectEqual(640, F[2].re);
    try std.testing.expectEqual(658, F[2].im);
    try std.testing.expectEqual(700, F[3].re);
    try std.testing.expectEqual(724, F[3].im);
    try std.testing.expectEqual(760, F[4].re);
    try std.testing.expectEqual(790, F[4].im);
    try std.testing.expectEqual(592, F[5].re);
    try std.testing.expectEqual(628, F[5].im);
    try std.testing.expectEqual(664, F[6].re);
    try std.testing.expectEqual(706, F[6].im);
    try std.testing.expectEqual(736, F[7].re);
    try std.testing.expectEqual(784, F[7].im);
    try std.testing.expectEqual(808, F[8].re);
    try std.testing.expectEqual(862, F[8].im);
    try std.testing.expectEqual(880, F[9].re);
    try std.testing.expectEqual(940, F[9].im);
    try std.testing.expectEqual(664, F[10].re);
    try std.testing.expectEqual(730, F[10].im);
    try std.testing.expectEqual(748, F[11].re);
    try std.testing.expectEqual(820, F[11].im);
    try std.testing.expectEqual(832, F[12].re);
    try std.testing.expectEqual(910, F[12].im);
    try std.testing.expectEqual(916, F[13].re);
    try std.testing.expectEqual(1000, F[13].im);
    try std.testing.expectEqual(1000, F[14].re);
    try std.testing.expectEqual(1090, F[14].im);
    try std.testing.expectEqual(736, F[15].re);
    try std.testing.expectEqual(832, F[15].im);
    try std.testing.expectEqual(832, F[16].re);
    try std.testing.expectEqual(934, F[16].im);
    try std.testing.expectEqual(928, F[17].re);
    try std.testing.expectEqual(1036, F[17].im);
    try std.testing.expectEqual(1024, F[18].re);
    try std.testing.expectEqual(1138, F[18].im);
    try std.testing.expectEqual(1120, F[19].re);
    try std.testing.expectEqual(1240, F[19].im);

    gemm(.col_major, .conj_trans, .no_trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(38, F[0].re);
    try std.testing.expectEqual(3194, F[0].im);
    try std.testing.expectEqual(92, F[1].re);
    try std.testing.expectEqual(3644, F[1].im);
    try std.testing.expectEqual(146, F[2].re);
    try std.testing.expectEqual(4094, F[2].im);
    try std.testing.expectEqual(200, F[3].re);
    try std.testing.expectEqual(4544, F[3].im);
    try std.testing.expectEqual(38, F[4].re);
    try std.testing.expectEqual(4778, F[4].im);
    try std.testing.expectEqual(200, F[5].re);
    try std.testing.expectEqual(3968, F[5].im);
    try std.testing.expectEqual(362, F[6].re);
    try std.testing.expectEqual(4598, F[6].im);
    try std.testing.expectEqual(524, F[7].re);
    try std.testing.expectEqual(5228, F[7].im);
    try std.testing.expectEqual(38, F[8].re);
    try std.testing.expectEqual(5210, F[8].im);
    try std.testing.expectEqual(308, F[9].re);
    try std.testing.expectEqual(5948, F[9].im);
    try std.testing.expectEqual(578, F[10].re);
    try std.testing.expectEqual(4958, F[10].im);
    try std.testing.expectEqual(848, F[11].re);
    try std.testing.expectEqual(5768, F[11].im);
    try std.testing.expectEqual(38, F[12].re);
    try std.testing.expectEqual(5498, F[12].im);
    try std.testing.expectEqual(416, F[13].re);
    try std.testing.expectEqual(6416, F[13].im);
    try std.testing.expectEqual(794, F[14].re);
    try std.testing.expectEqual(7334, F[14].im);
    try std.testing.expectEqual(1172, F[15].re);
    try std.testing.expectEqual(6164, F[15].im);
    try std.testing.expectEqual(38, F[16].re);
    try std.testing.expectEqual(5642, F[16].im);
    try std.testing.expectEqual(524, F[17].re);
    try std.testing.expectEqual(6740, F[17].im);
    try std.testing.expectEqual(1010, F[18].re);
    try std.testing.expectEqual(7838, F[18].im);
    try std.testing.expectEqual(1496, F[19].re);
    try std.testing.expectEqual(8936, F[19].im);

    gemm(.row_major, .conj_trans, .conj_no_trans, m, n, k, gamma, D.ptr, m, E.ptr, n, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-8948, F[0].re);
    try std.testing.expectEqual(9176, F[0].im);
    try std.testing.expectEqual(-10076, F[1].re);
    try std.testing.expectEqual(10628, F[1].im);
    try std.testing.expectEqual(-11204, F[2].re);
    try std.testing.expectEqual(12080, F[2].im);
    try std.testing.expectEqual(-12332, F[3].re);
    try std.testing.expectEqual(13532, F[3].im);
    try std.testing.expectEqual(-13460, F[4].re);
    try std.testing.expectEqual(13688, F[4].im);
    try std.testing.expectEqual(-10712, F[5].re);
    try std.testing.expectEqual(11912, F[5].im);
    try std.testing.expectEqual(-12044, F[6].re);
    try std.testing.expectEqual(14216, F[6].im);
    try std.testing.expectEqual(-13376, F[7].re);
    try std.testing.expectEqual(16520, F[7].im);
    try std.testing.expectEqual(-14708, F[8].re);
    try std.testing.expectEqual(14936, F[8].im);
    try std.testing.expectEqual(-16040, F[9].re);
    try std.testing.expectEqual(17888, F[9].im);
    try std.testing.expectEqual(-12476, F[10].re);
    try std.testing.expectEqual(15944, F[10].im);
    try std.testing.expectEqual(-14012, F[11].re);
    try std.testing.expectEqual(19100, F[11].im);
    try std.testing.expectEqual(-15548, F[12].re);
    try std.testing.expectEqual(15776, F[12].im);
    try std.testing.expectEqual(-17084, F[13].re);
    try std.testing.expectEqual(19580, F[13].im);
    try std.testing.expectEqual(-18620, F[14].re);
    try std.testing.expectEqual(23384, F[14].im);
    try std.testing.expectEqual(-14240, F[15].re);
    try std.testing.expectEqual(21272, F[15].im);
    try std.testing.expectEqual(-15980, F[16].re);
    try std.testing.expectEqual(16208, F[16].im);
    try std.testing.expectEqual(-17720, F[17].re);
    try std.testing.expectEqual(20864, F[17].im);
    try std.testing.expectEqual(-19460, F[18].re);
    try std.testing.expectEqual(25520, F[18].im);
    try std.testing.expectEqual(-21200, F[19].re);
    try std.testing.expectEqual(30176, F[19].im);

    gemm(.col_major, .conj_trans, .conj_no_trans, m, n, k, gamma, D.ptr, k, E.ptr, k, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-54316, F[0].re);
    try std.testing.expectEqual(628, F[0].im);
    try std.testing.expectEqual(-61984, F[1].re);
    try std.testing.expectEqual(1528, F[1].im);
    try std.testing.expectEqual(-69652, F[2].re);
    try std.testing.expectEqual(2428, F[2].im);
    try std.testing.expectEqual(-77320, F[3].re);
    try std.testing.expectEqual(3328, F[3].im);
    try std.testing.expectEqual(-81316, F[4].re);
    try std.testing.expectEqual(556, F[4].im);
    try std.testing.expectEqual(-67564, F[5].re);
    try std.testing.expectEqual(3292, F[5].im);
    try std.testing.expectEqual(-78292, F[6].re);
    try std.testing.expectEqual(6028, F[6].im);
    try std.testing.expectEqual(-89020, F[7].re);
    try std.testing.expectEqual(8764, F[7].im);
    try std.testing.expectEqual(-88732, F[8].re);
    try std.testing.expectEqual(484, F[8].im);
    try std.testing.expectEqual(-101296, F[9].re);
    try std.testing.expectEqual(5056, F[9].im);
    try std.testing.expectEqual(-84484, F[10].re);
    try std.testing.expectEqual(9628, F[10].im);
    try std.testing.expectEqual(-98272, F[11].re);
    try std.testing.expectEqual(14200, F[11].im);
    try std.testing.expectEqual(-93700, F[12].re);
    try std.testing.expectEqual(412, F[12].im);
    try std.testing.expectEqual(-109324, F[13].re);
    try std.testing.expectEqual(6820, F[13].im);
    try std.testing.expectEqual(-124948, F[14].re);
    try std.testing.expectEqual(13228, F[14].im);
    try std.testing.expectEqual(-105076, F[15].re);
    try std.testing.expectEqual(19636, F[15].im);
    try std.testing.expectEqual(-96220, F[16].re);
    try std.testing.expectEqual(340, F[16].im);
    try std.testing.expectEqual(-114904, F[17].re);
    try std.testing.expectEqual(8584, F[17].im);
    try std.testing.expectEqual(-133588, F[18].re);
    try std.testing.expectEqual(16828, F[18].im);
    try std.testing.expectEqual(-152272, F[19].re);
    try std.testing.expectEqual(25072, F[19].im);

    gemm(.row_major, .conj_trans, .trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-164680, F[0].re);
    try std.testing.expectEqual(-160912, F[0].im);
    try std.testing.expectEqual(-190204, F[1].re);
    try std.testing.expectEqual(-181036, F[1].im);
    try std.testing.expectEqual(-215728, F[2].re);
    try std.testing.expectEqual(-201160, F[2].im);
    try std.testing.expectEqual(-241252, F[3].re);
    try std.testing.expectEqual(-221284, F[3].im);
    try std.testing.expectEqual(-244744, F[4].re);
    try std.testing.expectEqual(-241408, F[4].im);
    try std.testing.expectEqual(-212392, F[5].re);
    try std.testing.expectEqual(-192640, F[5].im);
    try std.testing.expectEqual(-252568, F[6].re);
    try std.testing.expectEqual(-216400, F[6].im);
    try std.testing.expectEqual(-292744, F[7].re);
    try std.testing.expectEqual(-240160, F[7].im);
    try std.testing.expectEqual(-266824, F[8].re);
    try std.testing.expectEqual(-263920, F[8].im);
    try std.testing.expectEqual(-318016, F[9].re);
    try std.testing.expectEqual(-287680, F[9].im);
    try std.testing.expectEqual(-282136, F[10].re);
    try std.testing.expectEqual(-224368, F[10].im);
    try std.testing.expectEqual(-336964, F[11].re);
    try std.testing.expectEqual(-251764, F[11].im);
    try std.testing.expectEqual(-281632, F[12].re);
    try std.testing.expectEqual(-279160, F[12].im);
    try std.testing.expectEqual(-347476, F[13].re);
    try std.testing.expectEqual(-306556, F[13].im);
    try std.testing.expectEqual(-413320, F[14].re);
    try std.testing.expectEqual(-333952, F[14].im);
    try std.testing.expectEqual(-373912, F[15].re);
    try std.testing.expectEqual(-256096, F[15].im);
    try std.testing.expectEqual(-289168, F[16].re);
    try std.testing.expectEqual(-287128, F[16].im);
    try std.testing.expectEqual(-369664, F[17].re);
    try std.testing.expectEqual(-318160, F[17].im);
    try std.testing.expectEqual(-450160, F[18].re);
    try std.testing.expectEqual(-349192, F[18].im);
    try std.testing.expectEqual(-530656, F[19].re);
    try std.testing.expectEqual(-380224, F[19].im);

    gemm(.col_major, .conj_trans, .trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-11120, F[0].re);
    try std.testing.expectEqual(-976592, F[0].im);
    try std.testing.expectEqual(-27104, F[1].re);
    try std.testing.expectEqual(-1113320, F[1].im);
    try std.testing.expectEqual(-43088, F[2].re);
    try std.testing.expectEqual(-1250048, F[2].im);
    try std.testing.expectEqual(-59072, F[3].re);
    try std.testing.expectEqual(-1386776, F[3].im);
    try std.testing.expectEqual(-9800, F[4].re);
    try std.testing.expectEqual(-1458248, F[4].im);
    try std.testing.expectEqual(-58796, F[5].re);
    try std.testing.expectEqual(-1214636, F[5].im);
    try std.testing.expectEqual(-107792, F[6].re);
    try std.testing.expectEqual(-1406192, F[6].im);
    try std.testing.expectEqual(-156788, F[7].re);
    try std.testing.expectEqual(-1597748, F[7].im);
    try std.testing.expectEqual(-8480, F[8].re);
    try std.testing.expectEqual(-1592000, F[8].im);
    try std.testing.expectEqual(-90488, F[9].re);
    try std.testing.expectEqual(-1816568, F[9].im);
    try std.testing.expectEqual(-172496, F[10].re);
    try std.testing.expectEqual(-1518704, F[10].im);
    try std.testing.expectEqual(-254504, F[11].re);
    try std.testing.expectEqual(-1765088, F[11].im);
    try std.testing.expectEqual(-7160, F[12].re);
    try std.testing.expectEqual(-1682120, F[12].im);
    try std.testing.expectEqual(-122180, F[13].re);
    try std.testing.expectEqual(-1961516, F[13].im);
    try std.testing.expectEqual(-237200, F[14].re);
    try std.testing.expectEqual(-2240912, F[14].im);
    try std.testing.expectEqual(-352220, F[15].re);
    try std.testing.expectEqual(-1888796, F[15].im);
    try std.testing.expectEqual(-5840, F[16].re);
    try std.testing.expectEqual(-1728608, F[16].im);
    try std.testing.expectEqual(-153872, F[17].re);
    try std.testing.expectEqual(-2062832, F[17].im);
    try std.testing.expectEqual(-301904, F[18].re);
    try std.testing.expectEqual(-2397056, F[18].im);
    try std.testing.expectEqual(-449936, F[19].re);
    try std.testing.expectEqual(-2731280, F[19].im);

    gemm(.row_major, .conj_trans, .conj_trans, m, n, k, gamma, D.ptr, m, E.ptr, k, delta, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(2896568, F[0].re);
    try std.testing.expectEqual(-2963288, F[0].im);
    try std.testing.expectEqual(3258980, F[1].re);
    try std.testing.expectEqual(-3421604, F[1].im);
    try std.testing.expectEqual(3621392, F[2].re);
    try std.testing.expectEqual(-3879920, F[2].im);
    try std.testing.expectEqual(3983804, F[3].re);
    try std.testing.expectEqual(-4338236, F[3].im);
    try std.testing.expectEqual(4346216, F[4].re);
    try std.testing.expectEqual(-4405016, F[4].im);
    try std.testing.expectEqual(3467696, F[5].re);
    try std.testing.expectEqual(-3820472, F[5].im);
    try std.testing.expectEqual(3895592, F[6].re);
    try std.testing.expectEqual(-4542344, F[6].im);
    try std.testing.expectEqual(4323488, F[7].re);
    try std.testing.expectEqual(-5264216, F[7].im);
    try std.testing.expectEqual(4751384, F[8].re);
    try std.testing.expectEqual(-4802264, F[8].im);
    try std.testing.expectEqual(5179280, F[9].re);
    try std.testing.expectEqual(-5722208, F[9].im);
    try std.testing.expectEqual(4038824, F[10].re);
    try std.testing.expectEqual(-5073800, F[10].im);
    try std.testing.expectEqual(4532204, F[11].re);
    try std.testing.expectEqual(-6059228, F[11].im);
    try std.testing.expectEqual(5025584, F[12].re);
    try std.testing.expectEqual(-5068544, F[12].im);
    try std.testing.expectEqual(5518964, F[13].re);
    try std.testing.expectEqual(-6252044, F[13].im);
    try std.testing.expectEqual(6012344, F[14].re);
    try std.testing.expectEqual(-7435544, F[14].im);
    try std.testing.expectEqual(4609952, F[15].re);
    try std.testing.expectEqual(-6723272, F[15].im);
    try std.testing.expectEqual(5168816, F[16].re);
    try std.testing.expectEqual(-5203856, F[16].im);
    try std.testing.expectEqual(5727680, F[17].re);
    try std.testing.expectEqual(-6650912, F[17].im);
    try std.testing.expectEqual(6286544, F[18].re);
    try std.testing.expectEqual(-8097968, F[18].im);
    try std.testing.expectEqual(6845408, F[19].re);
    try std.testing.expectEqual(-9545024, F[19].im);

    gemm(.col_major, .conj_trans, .conj_trans, m, n, k, gamma, D.ptr, k, E.ptr, n, delta, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(17579752, F[0].re);
    try std.testing.expectEqual(-200344, F[0].im);
    try std.testing.expectEqual(20042152, F[1].re);
    try std.testing.expectEqual(-488272, F[1].im);
    try std.testing.expectEqual(22504552, F[2].re);
    try std.testing.expectEqual(-776200, F[2].im);
    try std.testing.expectEqual(24966952, F[3].re);
    try std.testing.expectEqual(-1064128, F[3].im);
    try std.testing.expectEqual(26253904, F[4].re);
    try std.testing.expectEqual(-176608, F[4].im);
    try std.testing.expectEqual(21864964, F[5].re);
    try std.testing.expectEqual(-1058788, F[5].im);
    try std.testing.expectEqual(25314520, F[6].re);
    try std.testing.expectEqual(-1940968, F[6].im);
    try std.testing.expectEqual(28764076, F[7].re);
    try std.testing.expectEqual(-2823148, F[7].im);
    try std.testing.expectEqual(28661176, F[8].re);
    try std.testing.expectEqual(-152872, F[8].im);
    try std.testing.expectEqual(32704984, F[9].re);
    try std.testing.expectEqual(-1629304, F[9].im);
    try std.testing.expectEqual(27338680, F[10].re);
    try std.testing.expectEqual(-3105736, F[10].im);
    try std.testing.expectEqual(31775392, F[11].re);
    try std.testing.expectEqual(-4582168, F[11].im);
    try std.testing.expectEqual(30282640, F[12].re);
    try std.testing.expectEqual(-129136, F[12].im);
    try std.testing.expectEqual(35313604, F[13].re);
    try std.testing.expectEqual(-2199820, F[13].im);
    try std.testing.expectEqual(40344568, F[14].re);
    try std.testing.expectEqual(-4270504, F[14].im);
    try std.testing.expectEqual(34000900, F[15].re);
    try std.testing.expectEqual(-6341188, F[15].im);
    try std.testing.expectEqual(31118296, F[16].re);
    try std.testing.expectEqual(-105400, F[16].im);
    try std.testing.expectEqual(37136416, F[17].re);
    try std.testing.expectEqual(-2770336, F[17].im);
    try std.testing.expectEqual(43154536, F[18].re);
    try std.testing.expectEqual(-5435272, F[18].im);
    try std.testing.expectEqual(49172656, F[19].re);
    try std.testing.expectEqual(-8100208, F[19].im);
}
