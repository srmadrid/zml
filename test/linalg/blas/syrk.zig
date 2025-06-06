const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const syrk = zml.linalg.blas.syrk;

test syrk {
    const a = std.testing.allocator;

    const n = 5;
    const k = 3;
    const alpha = 2;
    const beta = 3;

    const A = try a.alloc(f64, n * k);
    defer a.free(A);
    const B = try a.alloc(f64, n * n);
    defer a.free(B);

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

    syrk(f64, .RowMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(31, B[0]);
    try std.testing.expectEqual(70, B[1]);
    try std.testing.expectEqual(109, B[2]);
    try std.testing.expectEqual(148, B[3]);
    try std.testing.expectEqual(187, B[4]);
    try std.testing.expectEqual(6, B[5]);
    try std.testing.expectEqual(175, B[6]);
    try std.testing.expectEqual(268, B[7]);
    try std.testing.expectEqual(361, B[8]);
    try std.testing.expectEqual(454, B[9]);
    try std.testing.expectEqual(11, B[10]);
    try std.testing.expectEqual(12, B[11]);
    try std.testing.expectEqual(427, B[12]);
    try std.testing.expectEqual(574, B[13]);
    try std.testing.expectEqual(721, B[14]);
    try std.testing.expectEqual(16, B[15]);
    try std.testing.expectEqual(17, B[16]);
    try std.testing.expectEqual(18, B[17]);
    try std.testing.expectEqual(787, B[18]);
    try std.testing.expectEqual(988, B[19]);
    try std.testing.expectEqual(21, B[20]);
    try std.testing.expectEqual(22, B[21]);
    try std.testing.expectEqual(23, B[22]);
    try std.testing.expectEqual(24, B[23]);
    try std.testing.expectEqual(1255, B[24]);

    syrk(f64, .ColumnMajor, .Upper, .NoTrans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(409, B[0]);
    try std.testing.expectEqual(70, B[1]);
    try std.testing.expectEqual(109, B[2]);
    try std.testing.expectEqual(148, B[3]);
    try std.testing.expectEqual(187, B[4]);
    try std.testing.expectEqual(370, B[5]);
    try std.testing.expectEqual(919, B[6]);
    try std.testing.expectEqual(268, B[7]);
    try std.testing.expectEqual(361, B[8]);
    try std.testing.expectEqual(454, B[9]);
    try std.testing.expectEqual(421, B[10]);
    try std.testing.expectEqual(472, B[11]);
    try std.testing.expectEqual(1765, B[12]);
    try std.testing.expectEqual(574, B[13]);
    try std.testing.expectEqual(721, B[14]);
    try std.testing.expectEqual(472, B[15]);
    try std.testing.expectEqual(529, B[16]);
    try std.testing.expectEqual(586, B[17]);
    try std.testing.expectEqual(2947, B[18]);
    try std.testing.expectEqual(988, B[19]);
    try std.testing.expectEqual(523, B[20]);
    try std.testing.expectEqual(586, B[21]);
    try std.testing.expectEqual(649, B[22]);
    try std.testing.expectEqual(712, B[23]);
    try std.testing.expectEqual(4465, B[24]);

    syrk(f64, .RowMajor, .Upper, .Trans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(1543, B[0]);
    try std.testing.expectEqual(562, B[1]);
    try std.testing.expectEqual(715, B[2]);
    try std.testing.expectEqual(868, B[3]);
    try std.testing.expectEqual(1021, B[4]);
    try std.testing.expectEqual(370, B[5]);
    try std.testing.expectEqual(3151, B[6]);
    try std.testing.expectEqual(1240, B[7]);
    try std.testing.expectEqual(1561, B[8]);
    try std.testing.expectEqual(1882, B[9]);
    try std.testing.expectEqual(421, B[10]);
    try std.testing.expectEqual(472, B[11]);
    try std.testing.expectEqual(5779, B[12]);
    try std.testing.expectEqual(2254, B[13]);
    try std.testing.expectEqual(2743, B[14]);
    try std.testing.expectEqual(472, B[15]);
    try std.testing.expectEqual(529, B[16]);
    try std.testing.expectEqual(586, B[17]);
    try std.testing.expectEqual(9427, B[18]);
    try std.testing.expectEqual(3604, B[19]);
    try std.testing.expectEqual(523, B[20]);
    try std.testing.expectEqual(586, B[21]);
    try std.testing.expectEqual(649, B[22]);
    try std.testing.expectEqual(712, B[23]);
    try std.testing.expectEqual(14095, B[24]);

    syrk(f64, .ColumnMajor, .Upper, .Trans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(4657, B[0]);
    try std.testing.expectEqual(562, B[1]);
    try std.testing.expectEqual(715, B[2]);
    try std.testing.expectEqual(868, B[3]);
    try std.testing.expectEqual(1021, B[4]);
    try std.testing.expectEqual(1174, B[5]);
    try std.testing.expectEqual(9607, B[6]);
    try std.testing.expectEqual(1240, B[7]);
    try std.testing.expectEqual(1561, B[8]);
    try std.testing.expectEqual(1882, B[9]);
    try std.testing.expectEqual(1363, B[10]);
    try std.testing.expectEqual(1660, B[11]);
    try std.testing.expectEqual(17725, B[12]);
    try std.testing.expectEqual(2254, B[13]);
    try std.testing.expectEqual(2743, B[14]);
    try std.testing.expectEqual(1552, B[15]);
    try std.testing.expectEqual(1921, B[16]);
    try std.testing.expectEqual(2290, B[17]);
    try std.testing.expectEqual(29011, B[18]);
    try std.testing.expectEqual(3604, B[19]);
    try std.testing.expectEqual(1741, B[20]);
    try std.testing.expectEqual(2182, B[21]);
    try std.testing.expectEqual(2623, B[22]);
    try std.testing.expectEqual(3064, B[23]);
    try std.testing.expectEqual(43465, B[24]);

    syrk(f64, .RowMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(13999, B[0]);
    try std.testing.expectEqual(562, B[1]);
    try std.testing.expectEqual(715, B[2]);
    try std.testing.expectEqual(868, B[3]);
    try std.testing.expectEqual(1021, B[4]);
    try std.testing.expectEqual(3586, B[5]);
    try std.testing.expectEqual(28975, B[6]);
    try std.testing.expectEqual(1240, B[7]);
    try std.testing.expectEqual(1561, B[8]);
    try std.testing.expectEqual(1882, B[9]);
    try std.testing.expectEqual(4189, B[10]);
    try std.testing.expectEqual(5224, B[11]);
    try std.testing.expectEqual(53563, B[12]);
    try std.testing.expectEqual(2254, B[13]);
    try std.testing.expectEqual(2743, B[14]);
    try std.testing.expectEqual(4792, B[15]);
    try std.testing.expectEqual(6097, B[16]);
    try std.testing.expectEqual(7402, B[17]);
    try std.testing.expectEqual(87763, B[18]);
    try std.testing.expectEqual(3604, B[19]);
    try std.testing.expectEqual(5395, B[20]);
    try std.testing.expectEqual(6970, B[21]);
    try std.testing.expectEqual(8545, B[22]);
    try std.testing.expectEqual(10120, B[23]);
    try std.testing.expectEqual(131575, B[24]);

    syrk(f64, .ColumnMajor, .Lower, .NoTrans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(42313, B[0]);
    try std.testing.expectEqual(2038, B[1]);
    try std.testing.expectEqual(2533, B[2]);
    try std.testing.expectEqual(3028, B[3]);
    try std.testing.expectEqual(3523, B[4]);
    try std.testing.expectEqual(3586, B[5]);
    try std.testing.expectEqual(87319, B[6]);
    try std.testing.expectEqual(4156, B[7]);
    try std.testing.expectEqual(5161, B[8]);
    try std.testing.expectEqual(6166, B[9]);
    try std.testing.expectEqual(4189, B[10]);
    try std.testing.expectEqual(5224, B[11]);
    try std.testing.expectEqual(161173, B[12]);
    try std.testing.expectEqual(7294, B[13]);
    try std.testing.expectEqual(8809, B[14]);
    try std.testing.expectEqual(4792, B[15]);
    try std.testing.expectEqual(6097, B[16]);
    try std.testing.expectEqual(7402, B[17]);
    try std.testing.expectEqual(263875, B[18]);
    try std.testing.expectEqual(11452, B[19]);
    try std.testing.expectEqual(5395, B[20]);
    try std.testing.expectEqual(6970, B[21]);
    try std.testing.expectEqual(8545, B[22]);
    try std.testing.expectEqual(10120, B[23]);
    try std.testing.expectEqual(395425, B[24]);

    syrk(f64, .RowMajor, .Lower, .Trans, n, k, alpha, A.ptr, n, beta, B.ptr, n);

    try std.testing.expectEqual(127255, B[0]);
    try std.testing.expectEqual(2038, B[1]);
    try std.testing.expectEqual(2533, B[2]);
    try std.testing.expectEqual(3028, B[3]);
    try std.testing.expectEqual(3523, B[4]);
    try std.testing.expectEqual(11110, B[5]);
    try std.testing.expectEqual(262351, B[6]);
    try std.testing.expectEqual(4156, B[7]);
    try std.testing.expectEqual(5161, B[8]);
    try std.testing.expectEqual(6166, B[9]);
    try std.testing.expectEqual(12955, B[10]);
    try std.testing.expectEqual(16108, B[11]);
    try std.testing.expectEqual(484003, B[12]);
    try std.testing.expectEqual(7294, B[13]);
    try std.testing.expectEqual(8809, B[14]);
    try std.testing.expectEqual(14800, B[15]);
    try std.testing.expectEqual(18769, B[16]);
    try std.testing.expectEqual(22738, B[17]);
    try std.testing.expectEqual(792211, B[18]);
    try std.testing.expectEqual(11452, B[19]);
    try std.testing.expectEqual(16645, B[20]);
    try std.testing.expectEqual(21430, B[21]);
    try std.testing.expectEqual(26215, B[22]);
    try std.testing.expectEqual(31000, B[23]);
    try std.testing.expectEqual(1186975, B[24]);

    syrk(f64, .ColumnMajor, .Lower, .Trans, n, k, alpha, A.ptr, k, beta, B.ptr, n);

    try std.testing.expectEqual(381793, B[0]);
    try std.testing.expectEqual(6178, B[1]);
    try std.testing.expectEqual(7699, B[2]);
    try std.testing.expectEqual(9220, B[3]);
    try std.testing.expectEqual(10741, B[4]);
    try std.testing.expectEqual(11110, B[5]);
    try std.testing.expectEqual(787207, B[6]);
    try std.testing.expectEqual(12712, B[7]);
    try std.testing.expectEqual(15817, B[8]);
    try std.testing.expectEqual(18922, B[9]);
    try std.testing.expectEqual(12955, B[10]);
    try std.testing.expectEqual(16108, B[11]);
    try std.testing.expectEqual(1452397, B[12]);
    try std.testing.expectEqual(22414, B[13]);
    try std.testing.expectEqual(27103, B[14]);
    try std.testing.expectEqual(14800, B[15]);
    try std.testing.expectEqual(18769, B[16]);
    try std.testing.expectEqual(22738, B[17]);
    try std.testing.expectEqual(2377363, B[18]);
    try std.testing.expectEqual(35284, B[19]);
    try std.testing.expectEqual(16645, B[20]);
    try std.testing.expectEqual(21430, B[21]);
    try std.testing.expectEqual(26215, B[22]);
    try std.testing.expectEqual(31000, B[23]);
    try std.testing.expectEqual(3562105, B[24]);

    const gamma = cf64.init(2, 2);
    const delta = cf64.init(3, 3);

    const C = try a.alloc(cf64, n * k);
    defer a.free(C);
    const D = try a.alloc(cf64, n * n);
    defer a.free(D);

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
    });
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

    syrk(cf64, .RowMajor, .Upper, .NoTrans, n, k, gamma, C.ptr, k, delta, D.ptr, n);

    try std.testing.expectEqual(-56, D[0].re);
    try std.testing.expectEqual(62, D[0].im);
    try std.testing.expectEqual(-128, D[1].re);
    try std.testing.expectEqual(140, D[1].im);
    try std.testing.expectEqual(-200, D[2].re);
    try std.testing.expectEqual(218, D[2].im);
    try std.testing.expectEqual(-272, D[3].re);
    try std.testing.expectEqual(296, D[3].im);
    try std.testing.expectEqual(-344, D[4].re);
    try std.testing.expectEqual(374, D[4].im);
    try std.testing.expectEqual(6, D[5].re);
    try std.testing.expectEqual(6, D[5].im);
    try std.testing.expectEqual(-308, D[6].re);
    try std.testing.expectEqual(350, D[6].im);
    try std.testing.expectEqual(-488, D[7].re);
    try std.testing.expectEqual(536, D[7].im);
    try std.testing.expectEqual(-668, D[8].re);
    try std.testing.expectEqual(722, D[8].im);
    try std.testing.expectEqual(-848, D[9].re);
    try std.testing.expectEqual(908, D[9].im);
    try std.testing.expectEqual(11, D[10].re);
    try std.testing.expectEqual(11, D[10].im);
    try std.testing.expectEqual(12, D[11].re);
    try std.testing.expectEqual(12, D[11].im);
    try std.testing.expectEqual(-776, D[12].re);
    try std.testing.expectEqual(854, D[12].im);
    try std.testing.expectEqual(-1064, D[13].re);
    try std.testing.expectEqual(1148, D[13].im);
    try std.testing.expectEqual(-1352, D[14].re);
    try std.testing.expectEqual(1442, D[14].im);
    try std.testing.expectEqual(16, D[15].re);
    try std.testing.expectEqual(16, D[15].im);
    try std.testing.expectEqual(17, D[16].re);
    try std.testing.expectEqual(17, D[16].im);
    try std.testing.expectEqual(18, D[17].re);
    try std.testing.expectEqual(18, D[17].im);
    try std.testing.expectEqual(-1460, D[18].re);
    try std.testing.expectEqual(1574, D[18].im);
    try std.testing.expectEqual(-1856, D[19].re);
    try std.testing.expectEqual(1976, D[19].im);
    try std.testing.expectEqual(21, D[20].re);
    try std.testing.expectEqual(21, D[20].im);
    try std.testing.expectEqual(22, D[21].re);
    try std.testing.expectEqual(22, D[21].im);
    try std.testing.expectEqual(23, D[22].re);
    try std.testing.expectEqual(23, D[22].im);
    try std.testing.expectEqual(24, D[23].re);
    try std.testing.expectEqual(24, D[23].im);
    try std.testing.expectEqual(-2360, D[24].re);
    try std.testing.expectEqual(2510, D[24].im);

    syrk(cf64, .ColumnMajor, .Upper, .NoTrans, n, k, gamma, C.ptr, n, delta, D.ptr, n);

    try std.testing.expectEqual(-986, D[0].re);
    try std.testing.expectEqual(650, D[0].im);
    try std.testing.expectEqual(-128, D[1].re);
    try std.testing.expectEqual(140, D[1].im);
    try std.testing.expectEqual(-200, D[2].re);
    try std.testing.expectEqual(218, D[2].im);
    try std.testing.expectEqual(-272, D[3].re);
    try std.testing.expectEqual(296, D[3].im);
    try std.testing.expectEqual(-344, D[4].re);
    try std.testing.expectEqual(374, D[4].im);
    try std.testing.expectEqual(-704, D[5].re);
    try std.testing.expectEqual(740, D[5].im);
    try std.testing.expectEqual(-2762, D[6].re);
    try std.testing.expectEqual(914, D[6].im);
    try std.testing.expectEqual(-488, D[7].re);
    try std.testing.expectEqual(536, D[7].im);
    try std.testing.expectEqual(-668, D[8].re);
    try std.testing.expectEqual(722, D[8].im);
    try std.testing.expectEqual(-848, D[9].re);
    try std.testing.expectEqual(908, D[9].im);
    try std.testing.expectEqual(-776, D[10].re);
    try std.testing.expectEqual(842, D[10].im);
    try std.testing.expectEqual(-872, D[11].re);
    try std.testing.expectEqual(944, D[11].im);
    try std.testing.expectEqual(-5858, D[12].re);
    try std.testing.expectEqual(1202, D[12].im);
    try std.testing.expectEqual(-1064, D[13].re);
    try std.testing.expectEqual(1148, D[13].im);
    try std.testing.expectEqual(-1352, D[14].re);
    try std.testing.expectEqual(1442, D[14].im);
    try std.testing.expectEqual(-848, D[15].re);
    try std.testing.expectEqual(944, D[15].im);
    try std.testing.expectEqual(-956, D[16].re);
    try std.testing.expectEqual(1058, D[16].im);
    try std.testing.expectEqual(-1064, D[17].re);
    try std.testing.expectEqual(1172, D[17].im);
    try std.testing.expectEqual(-10274, D[18].re);
    try std.testing.expectEqual(1514, D[18].im);
    try std.testing.expectEqual(-1856, D[19].re);
    try std.testing.expectEqual(1976, D[19].im);
    try std.testing.expectEqual(-920, D[20].re);
    try std.testing.expectEqual(1046, D[20].im);
    try std.testing.expectEqual(-1040, D[21].re);
    try std.testing.expectEqual(1172, D[21].im);
    try std.testing.expectEqual(-1160, D[22].re);
    try std.testing.expectEqual(1298, D[22].im);
    try std.testing.expectEqual(-1280, D[23].re);
    try std.testing.expectEqual(1424, D[23].im);
    try std.testing.expectEqual(-16010, D[24].re);
    try std.testing.expectEqual(1850, D[24].im);

    syrk(cf64, .RowMajor, .Upper, .Trans, n, k, gamma, C.ptr, n, delta, D.ptr, n);

    try std.testing.expectEqual(-5540, D[0].re);
    try std.testing.expectEqual(-376, D[0].im);
    try std.testing.expectEqual(-1508, D[1].re);
    try std.testing.expectEqual(740, D[1].im);
    try std.testing.expectEqual(-2030, D[2].re);
    try std.testing.expectEqual(830, D[2].im);
    try std.testing.expectEqual(-2552, D[3].re);
    try std.testing.expectEqual(920, D[3].im);
    try std.testing.expectEqual(-3074, D[4].re);
    try std.testing.expectEqual(1010, D[4].im);
    try std.testing.expectEqual(-704, D[5].re);
    try std.testing.expectEqual(740, D[5].im);
    try std.testing.expectEqual(-11816, D[6].re);
    try std.testing.expectEqual(-4756, D[6].im);
    try std.testing.expectEqual(-3944, D[7].re);
    try std.testing.expectEqual(1016, D[7].im);
    try std.testing.expectEqual(-5126, D[8].re);
    try std.testing.expectEqual(1118, D[8].im);
    try std.testing.expectEqual(-6308, D[9].re);
    try std.testing.expectEqual(1220, D[9].im);
    try std.testing.expectEqual(-776, D[10].re);
    try std.testing.expectEqual(842, D[10].im);
    try std.testing.expectEqual(-872, D[11].re);
    try std.testing.expectEqual(944, D[11].im);
    try std.testing.expectEqual(-22148, D[12].re);
    try std.testing.expectEqual(-13000, D[12].im);
    try std.testing.expectEqual(-7700, D[13].re);
    try std.testing.expectEqual(1316, D[13].im);
    try std.testing.expectEqual(-9542, D[14].re);
    try std.testing.expectEqual(1430, D[14].im);
    try std.testing.expectEqual(-848, D[15].re);
    try std.testing.expectEqual(944, D[15].im);
    try std.testing.expectEqual(-956, D[16].re);
    try std.testing.expectEqual(1058, D[16].im);
    try std.testing.expectEqual(-1064, D[17].re);
    try std.testing.expectEqual(1172, D[17].im);
    try std.testing.expectEqual(-36536, D[18].re);
    try std.testing.expectEqual(-25108, D[18].im);
    try std.testing.expectEqual(-12776, D[19].re);
    try std.testing.expectEqual(1640, D[19].im);
    try std.testing.expectEqual(-920, D[20].re);
    try std.testing.expectEqual(1046, D[20].im);
    try std.testing.expectEqual(-1040, D[21].re);
    try std.testing.expectEqual(1172, D[21].im);
    try std.testing.expectEqual(-1160, D[22].re);
    try std.testing.expectEqual(1298, D[22].im);
    try std.testing.expectEqual(-1280, D[23].re);
    try std.testing.expectEqual(1424, D[23].im);
    try std.testing.expectEqual(-54980, D[24].re);
    try std.testing.expectEqual(-41080, D[24].im);

    syrk(cf64, .ColumnMajor, .Upper, .Trans, n, k, gamma, C.ptr, k, delta, D.ptr, n);

    try std.testing.expectEqual(-15548, D[0].re);
    try std.testing.expectEqual(-17692, D[0].im);
    try std.testing.expectEqual(-1508, D[1].re);
    try std.testing.expectEqual(740, D[1].im);
    try std.testing.expectEqual(-2030, D[2].re);
    try std.testing.expectEqual(830, D[2].im);
    try std.testing.expectEqual(-2552, D[3].re);
    try std.testing.expectEqual(920, D[3].im);
    try std.testing.expectEqual(-3074, D[4].re);
    try std.testing.expectEqual(1010, D[4].im);
    try std.testing.expectEqual(-4460, D[5].re);
    try std.testing.expectEqual(236, D[5].im);
    try std.testing.expectEqual(-21488, D[6].re);
    try std.testing.expectEqual(-49408, D[6].im);
    try std.testing.expectEqual(-3944, D[7].re);
    try std.testing.expectEqual(1016, D[7].im);
    try std.testing.expectEqual(-5126, D[8].re);
    try std.testing.expectEqual(1118, D[8].im);
    try std.testing.expectEqual(-6308, D[9].re);
    try std.testing.expectEqual(1220, D[9].im);
    try std.testing.expectEqual(-5054, D[10].re);
    try std.testing.expectEqual(398, D[10].im);
    try std.testing.expectEqual(-5936, D[11].re);
    try std.testing.expectEqual(704, D[11].im);
    try std.testing.expectEqual(-28220, D[12].re);
    try std.testing.expectEqual(-104668, D[12].im);
    try std.testing.expectEqual(-7700, D[13].re);
    try std.testing.expectEqual(1316, D[13].im);
    try std.testing.expectEqual(-9542, D[14].re);
    try std.testing.expectEqual(1430, D[14].im);
    try std.testing.expectEqual(-5648, D[15].re);
    try std.testing.expectEqual(560, D[15].im);
    try std.testing.expectEqual(-6710, D[16].re);
    try std.testing.expectEqual(974, D[16].im);
    try std.testing.expectEqual(-7772, D[17].re);
    try std.testing.expectEqual(1388, D[17].im);
    try std.testing.expectEqual(-35744, D[18].re);
    try std.testing.expectEqual(-183472, D[18].im);
    try std.testing.expectEqual(-12776, D[19].re);
    try std.testing.expectEqual(1640, D[19].im);
    try std.testing.expectEqual(-6242, D[20].re);
    try std.testing.expectEqual(722, D[20].im);
    try std.testing.expectEqual(-7484, D[21].re);
    try std.testing.expectEqual(1244, D[21].im);
    try std.testing.expectEqual(-8726, D[22].re);
    try std.testing.expectEqual(1766, D[22].im);
    try std.testing.expectEqual(-9968, D[23].re);
    try std.testing.expectEqual(2288, D[23].im);
    try std.testing.expectEqual(-44060, D[24].re);
    try std.testing.expectEqual(-285820, D[24].im);

    syrk(cf64, .RowMajor, .Lower, .NoTrans, n, k, gamma, C.ptr, k, delta, D.ptr, n);

    try std.testing.expectEqual(6376, D[0].re);
    try std.testing.expectEqual(-99664, D[0].im);
    try std.testing.expectEqual(-1508, D[1].re);
    try std.testing.expectEqual(740, D[1].im);
    try std.testing.expectEqual(-2030, D[2].re);
    try std.testing.expectEqual(830, D[2].im);
    try std.testing.expectEqual(-2552, D[3].re);
    try std.testing.expectEqual(920, D[3].im);
    try std.testing.expectEqual(-3074, D[4].re);
    try std.testing.expectEqual(1010, D[4].im);
    try std.testing.expectEqual(-14216, D[5].re);
    try std.testing.expectEqual(-12544, D[5].im);
    try std.testing.expectEqual(83452, D[6].re);
    try std.testing.expectEqual(-212380, D[6].im);
    try std.testing.expectEqual(-3944, D[7].re);
    try std.testing.expectEqual(1016, D[7].im);
    try std.testing.expectEqual(-5126, D[8].re);
    try std.testing.expectEqual(1118, D[8].im);
    try std.testing.expectEqual(-6308, D[9].re);
    try std.testing.expectEqual(1220, D[9].im);
    try std.testing.expectEqual(-16556, D[10].re);
    try std.testing.expectEqual(-13768, D[10].im);
    try std.testing.expectEqual(-20408, D[11].re);
    try std.testing.expectEqual(-15208, D[11].im);
    try std.testing.expectEqual(228568, D[12].re);
    try std.testing.expectEqual(-397888, D[12].im);
    try std.testing.expectEqual(-7700, D[13].re);
    try std.testing.expectEqual(1316, D[13].im);
    try std.testing.expectEqual(-9542, D[14].re);
    try std.testing.expectEqual(1430, D[14].im);
    try std.testing.expectEqual(-18896, D[15].re);
    try std.testing.expectEqual(-14992, D[15].im);
    try std.testing.expectEqual(-23720, D[16].re);
    try std.testing.expectEqual(-16540, D[16].im);
    try std.testing.expectEqual(-28544, D[17].re);
    try std.testing.expectEqual(-18088, D[17].im);
    try std.testing.expectEqual(441724, D[18].re);
    try std.testing.expectEqual(-656188, D[18].im);
    try std.testing.expectEqual(-12776, D[19].re);
    try std.testing.expectEqual(1640, D[19].im);
    try std.testing.expectEqual(-21236, D[20].re);
    try std.testing.expectEqual(-16216, D[20].im);
    try std.testing.expectEqual(-27032, D[21].re);
    try std.testing.expectEqual(-17872, D[21].im);
    try std.testing.expectEqual(-32828, D[22].re);
    try std.testing.expectEqual(-19528, D[22].im);
    try std.testing.expectEqual(-38624, D[23].re);
    try std.testing.expectEqual(-21184, D[23].im);
    try std.testing.expectEqual(722920, D[24].re);
    try std.testing.expectEqual(-987280, D[24].im);

    syrk(cf64, .ColumnMajor, .Lower, .NoTrans, n, k, gamma, C.ptr, n, delta, D.ptr, n);

    try std.testing.expectEqual(317488, D[0].re);
    try std.testing.expectEqual(-279232, D[0].im);
    try std.testing.expectEqual(-7448, D[1].re);
    try std.testing.expectEqual(-1600, D[1].im);
    try std.testing.expectEqual(-9356, D[2].re);
    try std.testing.expectEqual(-2824, D[2].im);
    try std.testing.expectEqual(-11264, D[3].re);
    try std.testing.expectEqual(-4048, D[3].im);
    try std.testing.expectEqual(-13172, D[4].re);
    try std.testing.expectEqual(-5272, D[4].im);
    try std.testing.expectEqual(-14216, D[5].re);
    try std.testing.expectEqual(-12544, D[5].im);
    try std.testing.expectEqual(886708, D[6].re);
    try std.testing.expectEqual(-385996, D[6].im);
    try std.testing.expectEqual(-15752, D[7].re);
    try std.testing.expectEqual(-7912, D[7].im);
    try std.testing.expectEqual(-19688, D[8].re);
    try std.testing.expectEqual(-11068, D[8].im);
    try std.testing.expectEqual(-23624, D[9].re);
    try std.testing.expectEqual(-14224, D[9].im);
    try std.testing.expectEqual(-16556, D[10].re);
    try std.testing.expectEqual(-13768, D[10].im);
    try std.testing.expectEqual(-20408, D[11].re);
    try std.testing.expectEqual(-15208, D[11].im);
    try std.testing.expectEqual(1878400, D[12].re);
    try std.testing.expectEqual(-506992, D[12].im);
    try std.testing.expectEqual(-28112, D[13].re);
    try std.testing.expectEqual(-18088, D[13].im);
    try std.testing.expectEqual(-34076, D[14].re);
    try std.testing.expectEqual(-23176, D[14].im);
    try std.testing.expectEqual(-18896, D[15].re);
    try std.testing.expectEqual(-14992, D[15].im);
    try std.testing.expectEqual(-23720, D[16].re);
    try std.testing.expectEqual(-16540, D[16].im);
    try std.testing.expectEqual(-28544, D[17].re);
    try std.testing.expectEqual(-18088, D[17].im);
    try std.testing.expectEqual(3292564, D[18].re);
    try std.testing.expectEqual(-642220, D[18].im);
    try std.testing.expectEqual(-44528, D[19].re);
    try std.testing.expectEqual(-32128, D[19].im);
    try std.testing.expectEqual(-21236, D[20].re);
    try std.testing.expectEqual(-16216, D[20].im);
    try std.testing.expectEqual(-27032, D[21].re);
    try std.testing.expectEqual(-17872, D[21].im);
    try std.testing.expectEqual(-32828, D[22].re);
    try std.testing.expectEqual(-19528, D[22].im);
    try std.testing.expectEqual(-38624, D[23].re);
    try std.testing.expectEqual(-21184, D[23].im);
    try std.testing.expectEqual(5129200, D[24].re);
    try std.testing.expectEqual(-791680, D[24].im);

    syrk(cf64, .RowMajor, .Lower, .Trans, n, k, gamma, C.ptr, n, delta, D.ptr, n);

    try std.testing.expectEqual(1789528, D[0].re);
    try std.testing.expectEqual(115400, D[0].im);
    try std.testing.expectEqual(-7448, D[1].re);
    try std.testing.expectEqual(-1600, D[1].im);
    try std.testing.expectEqual(-9356, D[2].re);
    try std.testing.expectEqual(-2824, D[2].im);
    try std.testing.expectEqual(-11264, D[3].re);
    try std.testing.expectEqual(-4048, D[3].im);
    try std.testing.expectEqual(-13172, D[4].re);
    try std.testing.expectEqual(-5272, D[4].im);
    try std.testing.expectEqual(-5720, D[5].re);
    try std.testing.expectEqual(-79576, D[5].im);
    try std.testing.expectEqual(3817324, D[6].re);
    try std.testing.expectEqual(1502924, D[6].im);
    try std.testing.expectEqual(-15752, D[7].re);
    try std.testing.expectEqual(-7912, D[7].im);
    try std.testing.expectEqual(-19688, D[8].re);
    try std.testing.expectEqual(-11068, D[8].im);
    try std.testing.expectEqual(-23624, D[9].re);
    try std.testing.expectEqual(-14224, D[9].im);
    try std.testing.expectEqual(-9140, D[10].re);
    try std.testing.expectEqual(-90196, D[10].im);
    try std.testing.expectEqual(-16472, D[11].re);
    try std.testing.expectEqual(-105976, D[11].im);
    try std.testing.expectEqual(7155208, D[12].re);
    try std.testing.expectEqual(4115192, D[12].im);
    try std.testing.expectEqual(-28112, D[13].re);
    try std.testing.expectEqual(-18088, D[13].im);
    try std.testing.expectEqual(-34076, D[14].re);
    try std.testing.expectEqual(-23176, D[14].im);
    try std.testing.expectEqual(-12560, D[15].re);
    try std.testing.expectEqual(-100816, D[15].im);
    try std.testing.expectEqual(-22496, D[16].re);
    try std.testing.expectEqual(-119824, D[16].im);
    try std.testing.expectEqual(-32432, D[17].re);
    try std.testing.expectEqual(-138832, D[17].im);
    try std.testing.expectEqual(11803180, D[18].re);
    try std.testing.expectEqual(7952204, D[18].im);
    try std.testing.expectEqual(-44528, D[19].re);
    try std.testing.expectEqual(-32128, D[19].im);
    try std.testing.expectEqual(-15980, D[20].re);
    try std.testing.expectEqual(-111436, D[20].im);
    try std.testing.expectEqual(-28520, D[21].re);
    try std.testing.expectEqual(-133672, D[21].im);
    try std.testing.expectEqual(-41060, D[22].re);
    try std.testing.expectEqual(-155908, D[22].im);
    try std.testing.expectEqual(-53600, D[23].re);
    try std.testing.expectEqual(-178144, D[23].im);
    try std.testing.expectEqual(17761240, D[24].re);
    try std.testing.expectEqual(13013960, D[24].im);

    syrk(cf64, .ColumnMajor, .Lower, .Trans, n, k, gamma, C.ptr, k, delta, D.ptr, n);

    try std.testing.expectEqual(5022328, D[0].re);
    try std.testing.expectEqual(5714840, D[0].im);
    try std.testing.expectEqual(-17672, D[1].re);
    try std.testing.expectEqual(-27016, D[1].im);
    try std.testing.expectEqual(-19796, D[2].re);
    try std.testing.expectEqual(-36340, D[2].im);
    try std.testing.expectEqual(-21920, D[3].re);
    try std.testing.expectEqual(-45664, D[3].im);
    try std.testing.expectEqual(-24044, D[4].re);
    try std.testing.expectEqual(-54988, D[4].im);
    try std.testing.expectEqual(-5720, D[5].re);
    try std.testing.expectEqual(-79576, D[5].im);
    try std.testing.expectEqual(6942892, D[6].re);
    try std.testing.expectEqual(15961052, D[6].im);
    try std.testing.expectEqual(-24008, D[7].re);
    try std.testing.expectEqual(-70504, D[7].im);
    try std.testing.expectEqual(-26528, D[8].re);
    try std.testing.expectEqual(-91600, D[8].im);
    try std.testing.expectEqual(-29048, D[9].re);
    try std.testing.expectEqual(-112696, D[9].im);
    try std.testing.expectEqual(-9140, D[10].re);
    try std.testing.expectEqual(-90196, D[10].im);
    try std.testing.expectEqual(-16472, D[11].re);
    try std.testing.expectEqual(-105976, D[11].im);
    try std.testing.expectEqual(9119272, D[12].re);
    try std.testing.expectEqual(33811976, D[12].im);
    try std.testing.expectEqual(-31136, D[13].re);
    try std.testing.expectEqual(-137536, D[13].im);
    try std.testing.expectEqual(-34052, D[14].re);
    try std.testing.expectEqual(-170404, D[14].im);
    try std.testing.expectEqual(-12560, D[15].re);
    try std.testing.expectEqual(-100816, D[15].im);
    try std.testing.expectEqual(-22496, D[16].re);
    try std.testing.expectEqual(-119824, D[16].im);
    try std.testing.expectEqual(-32432, D[17].re);
    try std.testing.expectEqual(-138832, D[17].im);
    try std.testing.expectEqual(11551468, D[18].re);
    try std.testing.expectEqual(59267612, D[18].im);
    try std.testing.expectEqual(-39056, D[19].re);
    try std.testing.expectEqual(-228112, D[19].im);
    try std.testing.expectEqual(-15980, D[20].re);
    try std.testing.expectEqual(-111436, D[20].im);
    try std.testing.expectEqual(-28520, D[21].re);
    try std.testing.expectEqual(-133672, D[21].im);
    try std.testing.expectEqual(-41060, D[22].re);
    try std.testing.expectEqual(-155908, D[22].im);
    try std.testing.expectEqual(-53600, D[23].re);
    try std.testing.expectEqual(-178144, D[23].im);
    try std.testing.expectEqual(14239480, D[24].re);
    try std.testing.expectEqual(92327960, D[24].im);
}
