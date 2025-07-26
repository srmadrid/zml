const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const trmm = zml.linalg.blas.trmm;

test trmm {
    const a = std.testing.allocator;

    const m = 4;
    const n = 5;
    const alpha = 2;

    const A = try a.alloc(f64, m * m);
    defer a.free(A);
    const B = try a.alloc(f64, n * n);
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
        13,
        14,
        15,
        16,
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

    trmm(.row_major, .left, .upper, .no_trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(220, C[0]);
    try std.testing.expectEqual(240, C[1]);
    try std.testing.expectEqual(260, C[2]);
    try std.testing.expectEqual(280, C[3]);
    try std.testing.expectEqual(300, C[4]);
    try std.testing.expectEqual(482, C[5]);
    try std.testing.expectEqual(524, C[6]);
    try std.testing.expectEqual(566, C[7]);
    try std.testing.expectEqual(608, C[8]);
    try std.testing.expectEqual(650, C[9]);
    try std.testing.expectEqual(626, C[10]);
    try std.testing.expectEqual(672, C[11]);
    try std.testing.expectEqual(718, C[12]);
    try std.testing.expectEqual(764, C[13]);
    try std.testing.expectEqual(810, C[14]);
    try std.testing.expectEqual(512, C[15]);
    try std.testing.expectEqual(544, C[16]);
    try std.testing.expectEqual(576, C[17]);
    try std.testing.expectEqual(608, C[18]);
    try std.testing.expectEqual(640, C[19]);

    trmm(.col_major, .left, .upper, .no_trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(14800, C[0]);
    try std.testing.expectEqual(15920, C[1]);
    try std.testing.expectEqual(14120, C[2]);
    try std.testing.expectEqual(8960, C[3]);
    try std.testing.expectEqual(29568, C[4]);
    try std.testing.expectEqual(32112, C[5]);
    try std.testing.expectEqual(28508, C[6]);
    try std.testing.expectEqual(18112, C[7]);
    try std.testing.expectEqual(36456, C[8]);
    try std.testing.expectEqual(39136, C[9]);
    try std.testing.expectEqual(33932, C[10]);
    try std.testing.expectEqual(21504, C[11]);
    try std.testing.expectEqual(36968, C[12]);
    try std.testing.expectEqual(39704, C[13]);
    try std.testing.expectEqual(33180, C[14]);
    try std.testing.expectEqual(16384, C[15]);
    try std.testing.expectEqual(34432, C[16]);
    try std.testing.expectEqual(36992, C[17]);
    try std.testing.expectEqual(32576, C[18]);
    try std.testing.expectEqual(20480, C[19]);

    trmm(.row_major, .left, .upper, .no_trans, .unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(492712, C[0]);
    try std.testing.expectEqual(550352, C[1]);
    try std.testing.expectEqual(618432, C[2]);
    try std.testing.expectEqual(662576, C[3]);
    try std.testing.expectEqual(578600, C[4]);
    try std.testing.expectEqual(801416, C[5]);
    try std.testing.expectEqual(908984, C[6]);
    try std.testing.expectEqual(1145648, C[7]);
    try std.testing.expectEqual(1149984, C[8]);
    try std.testing.expectEqual(870472, C[9]);
    try std.testing.expectEqual(461080, C[10]);
    try std.testing.expectEqual(869376, C[11]);
    try std.testing.expectEqual(961744, C[12]);
    try std.testing.expectEqual(861232, C[13]);
    try std.testing.expectEqual(557880, C[14]);
    try std.testing.expectEqual(32768, C[15]);
    try std.testing.expectEqual(68864, C[16]);
    try std.testing.expectEqual(73984, C[17]);
    try std.testing.expectEqual(65152, C[18]);
    try std.testing.expectEqual(40960, C[19]);

    trmm(.col_major, .left, .upper, .no_trans, .unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(34847696, C[0]);
    try std.testing.expectEqual(32021472, C[1]);
    try std.testing.expectEqual(21114144, C[2]);
    try std.testing.expectEqual(1325152, C[3]);
    try std.testing.expectEqual(55319920, C[4]);
    try std.testing.expectEqual(51860656, C[5]);
    try std.testing.expectEqual(36187408, C[6]);
    try std.testing.expectEqual(2291296, C[7]);
    try std.testing.expectEqual(41907904, C[8]);
    try std.testing.expectEqual(35305072, C[9]);
    try std.testing.expectEqual(27003440, C[10]);
    try std.testing.expectEqual(1738752, C[11]);
    try std.testing.expectEqual(21429616, C[12]);
    try std.testing.expectEqual(13797568, C[13]);
    try std.testing.expectEqual(2098800, C[14]);
    try std.testing.expectEqual(65536, C[15]);
    try std.testing.expectEqual(3115264, C[16]);
    try std.testing.expectEqual(2597888, C[17]);
    try std.testing.expectEqual(1359104, C[18]);
    try std.testing.expectEqual(81920, C[19]);

    trmm(.row_major, .left, .upper, .trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(69695392, C[0]);
    try std.testing.expectEqual(64042944, C[1]);
    try std.testing.expectEqual(42228288, C[2]);
    try std.testing.expectEqual(2650304, C[3]);
    try std.testing.expectEqual(110639840, C[4]);
    try std.testing.expectEqual(761718656, C[5]);
    try std.testing.expectEqual(562334784, C[6]);
    try std.testing.expectEqual(111952128, C[7]);
    try std.testing.expectEqual(508195456, C[8]);
    try std.testing.expectEqual(644940544, C[9]);
    try std.testing.expectEqual(1529211040, C[10]);
    try std.testing.expectEqual(737005088, C[11]);
    try std.testing.expectEqual(630214560, C[12]);
    try std.testing.expectEqual(898208064, C[13]);
    try std.testing.expectEqual(872364128, C[14]);
    try std.testing.expectEqual(1758731776, C[15]);
    try std.testing.expectEqual(976588800, C[16]);
    try std.testing.expectEqual(803017088, C[17]);
    try std.testing.expectEqual(1055760640, C[18]);
    try std.testing.expectEqual(1060433152, C[19]);

    trmm(.col_major, .left, .upper, .trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(139390784, C[0]);
    try std.testing.expectEqual(1465469248, C[1]);
    try std.testing.expectEqual(3464398272, C[2]);
    try std.testing.expectEqual(4956940992, C[3]);
    try std.testing.expectEqual(221279680, C[4]);
    try std.testing.expectEqual(10247022272, C[5]);
    try std.testing.expectEqual(29597255488, C[6]);
    try std.testing.expectEqual(44657269824, C[7]);
    try std.testing.expectEqual(1016390912, C[8]);
    try std.testing.expectEqual(12821241088, C[9]);
    try std.testing.expectEqual(55688971968, C[10]);
    try std.testing.expectEqual(100731911104, C[11]);
    try std.testing.expectEqual(1260429120, C[12]);
    try std.testing.expectEqual(17080642368, C[13]);
    try std.testing.expectEqual(48500034176, C[14]);
    try std.testing.expectEqual(123985745024, C[15]);
    try std.testing.expectEqual(1953177600, C[16]);
    try std.testing.expectEqual(19402093056, C[17]);
    try std.testing.expectEqual(56865674240, C[18]);
    try std.testing.expectEqual(113482467328, C[19]);

    trmm(.row_major, .left, .upper, .trans, .unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(278781568, C[0]);
    try std.testing.expectEqual(2930938496, C[1]);
    try std.testing.expectEqual(6928796544, C[2]);
    try std.testing.expectEqual(9913881984, C[3]);
    try std.testing.expectEqual(442559360, C[4]);
    try std.testing.expectEqual(21051607680, C[5]);
    try std.testing.expectEqual(65056387968, C[6]);
    try std.testing.expectEqual(103172132736, C[7]);
    try std.testing.expectEqual(21860545792, C[8]);
    try std.testing.expectEqual(26527600896, C[9]);
    try std.testing.expectEqual(255672600448, C[10]);
    try std.testing.expectEqual(624618214528, C[11]);
    try std.testing.expectEqual(648509025408, C[12]);
    try std.testing.expectEqual(78132403456, C[13]);
    try std.testing.expectEqual(277825121664, C[14]);
    try std.testing.expectEqual(1749574299904, C[15]);
    try std.testing.expectEqual(2906752063488, C[16]);
    try std.testing.expectEqual(811285988352, C[17]);
    try std.testing.expectEqual(579584547840, C[18]);
    try std.testing.expectEqual(1597875849728, C[19]);

    trmm(.col_major, .left, .upper, .trans, .unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(557563136, C[0]);
    try std.testing.expectEqual(8649692672, C[1]);
    try std.testing.expectEqual(77494431232, C[2]);
    try std.testing.expectEqual(317006258944, C[3]);
    try std.testing.expectEqual(885118720, C[4]);
    try std.testing.expectEqual(46528808960, C[5]);
    try std.testing.expectEqual(559110998016, C[6]);
    try std.testing.expectEqual(2758987462912, C[7]);
    try std.testing.expectEqual(43721091584, C[8]);
    try std.testing.expectEqual(271660659712, C[9]);
    try std.testing.expectEqual(1435387043072, C[10]);
    try std.testing.expectEqual(10230561458176, C[11]);
    try std.testing.expectEqual(1297018050816, C[12]);
    try std.testing.expectEqual(6641355060992, C[13]);
    try std.testing.expectEqual(13791460769792, C[14]);
    try std.testing.expectEqual(30882844207104, C[15]);
    try std.testing.expectEqual(5813504126976, C[16]);
    try std.testing.expectEqual(30690092611584, C[17]);
    try std.testing.expectEqual(69706426005504, C[18]);
    try std.testing.expectEqual(118874849459200, C[19]);

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

    trmm(.row_major, .left, .lower, .no_trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(2, C[0]);
    try std.testing.expectEqual(4, C[1]);
    try std.testing.expectEqual(6, C[2]);
    try std.testing.expectEqual(8, C[3]);
    try std.testing.expectEqual(10, C[4]);
    try std.testing.expectEqual(82, C[5]);
    try std.testing.expectEqual(104, C[6]);
    try std.testing.expectEqual(126, C[7]);
    try std.testing.expectEqual(148, C[8]);
    try std.testing.expectEqual(170, C[9]);
    try std.testing.expectEqual(380, C[10]);
    try std.testing.expectEqual(440, C[11]);
    try std.testing.expectEqual(500, C[12]);
    try std.testing.expectEqual(560, C[13]);
    try std.testing.expectEqual(620, C[14]);
    try std.testing.expectEqual(1036, C[15]);
    try std.testing.expectEqual(1152, C[16]);
    try std.testing.expectEqual(1268, C[17]);
    try std.testing.expectEqual(1384, C[18]);
    try std.testing.expectEqual(1500, C[19]);

    trmm(.col_major, .left, .lower, .no_trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(4, C[0]);
    try std.testing.expectEqual(56, C[1]);
    try std.testing.expectEqual(200, C[2]);
    try std.testing.expectEqual(480, C[3]);
    try std.testing.expectEqual(20, C[4]);
    try std.testing.expectEqual(1024, C[5]);
    try std.testing.expectEqual(3496, C[6]);
    try std.testing.expectEqual(7920, C[7]);
    try std.testing.expectEqual(296, C[8]);
    try std.testing.expectEqual(2632, C[9]);
    try std.testing.expectEqual(11628, C[10]);
    try std.testing.expectEqual(27104, C[11]);
    try std.testing.expectEqual(1000, C[12]);
    try std.testing.expectEqual(8720, C[13]);
    try std.testing.expectEqual(24480, C[14]);
    try std.testing.expectEqual(60992, C[15]);
    try std.testing.expectEqual(2304, C[16]);
    try std.testing.expectEqual(19824, C[17]);
    try std.testing.expectEqual(55112, C[18]);
    try std.testing.expectEqual(110720, C[19]);

    trmm(.row_major, .left, .lower, .no_trans, .unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(8, C[0]);
    try std.testing.expectEqual(112, C[1]);
    try std.testing.expectEqual(400, C[2]);
    try std.testing.expectEqual(960, C[3]);
    try std.testing.expectEqual(40, C[4]);
    try std.testing.expectEqual(2088, C[5]);
    try std.testing.expectEqual(7552, C[6]);
    try std.testing.expectEqual(17840, C[7]);
    try std.testing.expectEqual(5392, C[8]);
    try std.testing.expectEqual(5464, C[9]);
    try std.testing.expectEqual(43808, C[10]);
    try std.testing.expectEqual(125136, C[11]);
    try std.testing.expectEqual(164000, C[12]);
    try std.testing.expectEqual(32000, C[13]);
    try std.testing.expectEqual(101960, C[14]);
    try std.testing.expectEqual(499600, C[15]);
    try std.testing.expectEqual(917072, C[16]);
    try std.testing.expectEqual(296608, C[17]);
    try std.testing.expectEqual(392592, C[18]);
    try std.testing.expectEqual(1030056, C[19]);

    trmm(.col_major, .left, .lower, .no_trans, .unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(16, C[0]);
    try std.testing.expectEqual(256, C[1]);
    try std.testing.expectEqual(2416, C[2]);
    try std.testing.expectEqual(13376, C[3]);
    try std.testing.expectEqual(80, C[4]);
    try std.testing.expectEqual(4336, C[5]);
    try std.testing.expectEqual(44576, C[6]);
    try std.testing.expectEqual(250656, C[7]);
    try std.testing.expectEqual(10784, C[8]);
    try std.testing.expectEqual(32496, C[9]);
    try std.testing.expectEqual(196464, C[10]);
    try std.testing.expectEqual(1432224, C[11]);
    try std.testing.expectEqual(328000, C[12]);
    try std.testing.expectEqual(720000, C[13]);
    try std.testing.expectEqual(1635920, C[14]);
    try std.testing.expectEqual(5270240, C[15]);
    try std.testing.expectEqual(1834144, C[16]);
    try std.testing.expectEqual(4261504, C[17]);
    try std.testing.expectEqual(10440128, C[18]);
    try std.testing.expectEqual(23564624, C[19]);

    trmm(.row_major, .left, .lower, .trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(140605984, C[0]);
    try std.testing.expectEqual(73914048, C[1]);
    try std.testing.expectEqual(119214496, C[2]);
    try std.testing.expectEqual(284537920, C[3]);
    try std.testing.expectEqual(642451904, C[4]);
    try std.testing.expectEqual(151548032, C[5]);
    try std.testing.expectEqual(80535424, C[6]);
    try std.testing.expectEqual(128889984, C[7]);
    try std.testing.expectEqual(306852992, C[8]);
    try std.testing.expectEqual(692917824, C[9]);
    try std.testing.expectEqual(162429408, C[10]);
    try std.testing.expectEqual(86533248, C[11]);
    try std.testing.expectEqual(135061120, C[12]);
    try std.testing.expectEqual(329043840, C[13]);
    try std.testing.expectEqual(742928960, C[14]);
    try std.testing.expectEqual(168647680, C[15]);
    try std.testing.expectEqual(58692608, C[16]);
    try std.testing.expectEqual(136368128, C[17]);
    try std.testing.expectEqual(334084096, C[18]);
    try std.testing.expectEqual(754067968, C[19]);

    trmm(.col_major, .left, .lower, .trans, .non_unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(3568458496, C[0]);
    try std.testing.expectEqual(7108578240, C[1]);
    try std.testing.expectEqual(9451628992, C[2]);
    try std.testing.expectEqual(9105213440, C[3]);
    try std.testing.expectEqual(3405428352, C[4]);
    try std.testing.expectEqual(5008312064, C[5]);
    try std.testing.expectEqual(4865138944, C[6]);
    try std.testing.expectEqual(4124479488, C[7]);
    try std.testing.expectEqual(5052219712, C[8]);
    try std.testing.expectEqual(11973557568, C[9]);
    try std.testing.expectEqual(5650244928, C[10]);
    try std.testing.expectEqual(2769063936, C[11]);
    try std.testing.expectEqual(7393052800, C[12]);
    try std.testing.expectEqual(17047894400, C[13]);
    try std.testing.expectEqual(20391981440, C[14]);
    try std.testing.expectEqual(5396725760, C[15]);
    try std.testing.expectEqual(8699906048, C[16]);
    try std.testing.expectEqual(18378682368, C[17]);
    try std.testing.expectEqual(25447481344, C[18]);
    try std.testing.expectEqual(24130174976, C[19]);

    trmm(.row_major, .left, .lower, .trans, .unit, m, n, alpha, A.ptr, m, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(299239316096, C[0]);
    try std.testing.expectEqual(338909254016, C[1]);
    try std.testing.expectEqual(671068744832, C[2]);
    try std.testing.expectEqual(1037229238144, C[3]);
    try std.testing.expectEqual(1120986647680, C[4]);
    try std.testing.expectEqual(274129843968, C[5]);
    try std.testing.expectEqual(308708925952, C[6]);
    try std.testing.expectEqual(670713121280, C[7]);
    try std.testing.expectEqual(1063591805056, C[8]);
    try std.testing.expectEqual(1107431643264, C[9]);
    try std.testing.expectEqual(173202262656, C[10]);
    try std.testing.expectEqual(266535309312, C[11]);
    try std.testing.expectEqual(566146576640, C[12]);
    try std.testing.expectEqual(797520229120, C[13]);
    try std.testing.expectEqual(764689212160, C[14]);
    try std.testing.expectEqual(10793451520, C[15]);
    try std.testing.expectEqual(17399812096, C[16]);
    try std.testing.expectEqual(36757364736, C[17]);
    try std.testing.expectEqual(50894962688, C[18]);
    try std.testing.expectEqual(48260349952, C[19]);

    trmm(.col_major, .left, .lower, .trans, .unit, m, n, alpha, A.ptr, m, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(14278362022400, C[0]);
    try std.testing.expectEqual(26668448745984, C[1]);
    try std.testing.expectEqual(26235639205120, C[2]);
    try std.testing.expectEqual(2074458476288, C[3]);
    try std.testing.expectEqual(10556451197184, C[4]);
    try std.testing.expectEqual(15601594591744, C[5]);
    try std.testing.expectEqual(16714532762624, C[6]);
    try std.testing.expectEqual(1341426242560, C[7]);
    try std.testing.expectEqual(9728406233600, C[8]);
    try std.testing.expectEqual(8904259912704, C[9]);
    try std.testing.expectEqual(6743251948800, C[10]);
    try std.testing.expectEqual(533070618624, C[11]);
    try std.testing.expectEqual(8996856954880, C[12]);
    try std.testing.expectEqual(12473384652800, C[13]);
    try std.testing.expectEqual(1788421260800, C[14]);
    try std.testing.expectEqual(21586903040, C[15]);
    try std.testing.expectEqual(873281658880, C[16]);
    try std.testing.expectEqual(1558209806336, C[17]);
    try std.testing.expectEqual(1260038324224, C[18]);
    try std.testing.expectEqual(96520699904, C[19]);

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

    trmm(.row_major, .right, .upper, .no_trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(2, C[0]);
    try std.testing.expectEqual(32, C[1]);
    try std.testing.expectEqual(116, C[2]);
    try std.testing.expectEqual(280, C[3]);
    try std.testing.expectEqual(550, C[4]);
    try std.testing.expectEqual(12, C[5]);
    try std.testing.expectEqual(122, C[6]);
    try std.testing.expectEqual(356, C[7]);
    try std.testing.expectEqual(740, C[8]);
    try std.testing.expectEqual(1300, C[9]);
    try std.testing.expectEqual(22, C[10]);
    try std.testing.expectEqual(212, C[11]);
    try std.testing.expectEqual(596, C[12]);
    try std.testing.expectEqual(1200, C[13]);
    try std.testing.expectEqual(2050, C[14]);
    try std.testing.expectEqual(32, C[15]);
    try std.testing.expectEqual(302, C[16]);
    try std.testing.expectEqual(836, C[17]);
    try std.testing.expectEqual(1660, C[18]);
    try std.testing.expectEqual(2800, C[19]);

    trmm(.col_major, .right, .upper, .no_trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(4, C[0]);
    try std.testing.expectEqual(64, C[1]);
    try std.testing.expectEqual(232, C[2]);
    try std.testing.expectEqual(560, C[3]);
    try std.testing.expectEqual(7724, C[4]);
    try std.testing.expectEqual(552, C[5]);
    try std.testing.expectEqual(3100, C[6]);
    try std.testing.expectEqual(8344, C[7]);
    try std.testing.expectEqual(32484, C[8]);
    try std.testing.expectEqual(34792, C[9]);
    try std.testing.expectEqual(6052, C[10]);
    try std.testing.expectEqual(20216, C[11]);
    try std.testing.expectEqual(68052, C[12]);
    try std.testing.expectEqual(93832, C[13]);
    try std.testing.expectEqual(86552, C[14]);
    try std.testing.expectEqual(29912, C[15]);
    try std.testing.expectEqual(102032, C[16]);
    try std.testing.expectEqual(161072, C[17]);
    try std.testing.expectEqual(192652, C[18]);
    try std.testing.expectEqual(178712, C[19]);

    trmm(.row_major, .right, .upper, .no_trans, .unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(8, C[0]);
    try std.testing.expectEqual(144, C[1]);
    try std.testing.expectEqual(1512, C[2]);
    try std.testing.expectEqual(8800, C[3]);
    try std.testing.expectEqual(46128, C[4]);
    try std.testing.expectEqual(1104, C[5]);
    try std.testing.expectEqual(8408, C[6]);
    try std.testing.expectEqual(69600, C[7]);
    try std.testing.expectEqual(358816, C[8]);
    try std.testing.expectEqual(1686784, C[9]);
    try std.testing.expectEqual(12104, C[10]);
    try std.testing.expectEqual(64640, C[11]);
    try std.testing.expectEqual(495872, C[12]);
    try std.testing.expectEqual(2505424, C[13]);
    try std.testing.expectEqual(6432784, C[14]);
    try std.testing.expectEqual(59824, C[15]);
    try std.testing.expectEqual(323712, C[16]);
    try std.testing.expectEqual(2134128, C[17]);
    try std.testing.expectEqual(6971192, C[18]);
    try std.testing.expectEqual(15235424, C[19]);

    trmm(.col_major, .right, .upper, .no_trans, .unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(16, C[0]);
    try std.testing.expectEqual(288, C[1]);
    try std.testing.expectEqual(3024, C[2]);
    try std.testing.expectEqual(17600, C[3]);
    try std.testing.expectEqual(92352, C[4]);
    try std.testing.expectEqual(3936, C[5]);
    try std.testing.expectEqual(34960, C[6]);
    try std.testing.expectEqual(244800, C[7]);
    try std.testing.expectEqual(1824880, C[8]);
    try std.testing.expectEqual(3403232, C[9]);
    try std.testing.expectEqual(259264, C[10]);
    try std.testing.expectEqual(1993280, C[11]);
    try std.testing.expectEqual(15477728, C[12]);
    try std.testing.expectEqual(65777216, C[13]);
    try std.testing.expectEqual(13635568, C[14]);
    try std.testing.expectEqual(5094688, C[15]);
    try std.testing.expectEqual(42984784, C[16]);
    try std.testing.expectEqual(202175296, C[17]);
    try std.testing.expectEqual(323706256, C[18]);
    try std.testing.expectEqual(39747840, C[19]);

    trmm(.row_major, .right, .upper, .trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(1083648, C[0]);
    try std.testing.expectEqual(2216256, C[1]);
    try std.testing.expectEqual(3341984, C[2]);
    try std.testing.expectEqual(4362880, C[3]);
    try std.testing.expectEqual(4617600, C[4]);
    try std.testing.expectEqual(50247872, C[5]);
    try std.testing.expectEqual(105318720, C[6]);
    try std.testing.expectEqual(159558400, C[7]);
    try std.testing.expectEqual(205474720, C[8]);
    try std.testing.expectEqual(170161600, C[9]);
    try std.testing.expectEqual(763931424, C[10]);
    try std.testing.expectEqual(1732250816, C[11]);
    try std.testing.expectEqual(2653250016, C[12]);
    try std.testing.expectEqual(3044956928, C[13]);
    try std.testing.expectEqual(681778400, C[14]);
    try std.testing.expectEqual(4382308736, C[15]);
    try std.testing.expectEqual(10458261120, C[16]);
    try std.testing.expectEqual(15512768064, C[17]);
    try std.testing.expectEqual(13890751328, C[18]);
    try std.testing.expectEqual(1987392000, C[19]);

    trmm(.col_major, .right, .upper, .trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(528728989888, C[0]);
    try std.testing.expectEqual(753325842560, C[1]);
    try std.testing.expectEqual(623305464512, C[2]);
    try std.testing.expectEqual(263737288064, C[3]);
    try std.testing.expectEqual(555370029504, C[4]);
    try std.testing.expectEqual(790877678976, C[5]);
    try std.testing.expectEqual(654182340288, C[6]);
    try std.testing.expectEqual(280251582208, C[7]);
    try std.testing.expectEqual(581939354816, C[8]);
    try std.testing.expectEqual(827629981952, C[9]);
    try std.testing.expectEqual(683380800512, C[10]);
    try std.testing.expectEqual(294221667712, C[11]);
    try std.testing.expectEqual(602820034368, C[12]);
    try std.testing.expectEqual(860321230336, C[13]);
    try std.testing.expectEqual(692663642944, C[14]);
    try std.testing.expectEqual(261922547968, C[15]);
    try std.testing.expectEqual(522913056000, C[16]);
    try std.testing.expectEqual(775638403200, C[17]);
    try std.testing.expectEqual(694537566400, C[18]);
    try std.testing.expectEqual(99369600000, C[19]);

    trmm(.row_major, .right, .upper, .trans, .unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(15474192736640, C[0]);
    try std.testing.expectEqual(27334210892544, C[1]);
    try std.testing.expectEqual(25292355879936, C[2]);
    try std.testing.expectEqual(22742275756288, C[3]);
    try std.testing.expectEqual(1110740059008, C[4]);
    try std.testing.expectEqual(18811808870400, C[5]);
    try std.testing.expectEqual(32819898021632, C[6]);
    try std.testing.expectEqual(41683704557824, C[7]);
    try std.testing.expectEqual(34269077987712, C[8]);
    try std.testing.expectEqual(1655259963904, C[9]);
    try std.testing.expectEqual(19969774750208, C[10]);
    try std.testing.expectEqual(39572618890240, C[11]);
    try std.testing.expectEqual(46074543806464, C[12]);
    try std.testing.expectEqual(29427188178432, C[13]);
    try std.testing.expectEqual(1385327285888, C[14]);
    try std.testing.expectEqual(13819324270336, C[15]);
    try std.testing.expectEqual(27945108758400, C[16]);
    try std.testing.expectEqual(23979416665600, C[17]);
    try std.testing.expectEqual(5363859132800, C[18]);
    try std.testing.expectEqual(198739200000, C[19]);

    trmm(.col_major, .right, .upper, .trans, .unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(3446276951570688, C[0]);
    try std.testing.expectEqual(2265631369100800, C[1]);
    try std.testing.expectEqual(1153371089250048, C[2]);
    try std.testing.expectEqual(1866852044842496, C[3]);
    try std.testing.expectEqual(3620798626612480, C[4]);
    try std.testing.expectEqual(2132968588227584, C[5]);
    try std.testing.expectEqual(828025319611648, C[6]);
    try std.testing.expectEqual(1511711812472832, C[7]);
    try std.testing.expectEqual(3012696735894528, C[8]);
    try std.testing.expectEqual(2165742460968960, C[9]);
    try std.testing.expectEqual(336548851901184, C[10]);
    try std.testing.expectEqual(585782914712576, C[11]);
    try std.testing.expectEqual(1433514308016128, C[12]);
    try std.testing.expectEqual(1209866376305664, C[13]);
    try std.testing.expectEqual(260235892946176, C[14]);
    try std.testing.expectEqual(37178130140672, C[15]);
    try std.testing.expectEqual(55890217516800, C[16]);
    try std.testing.expectEqual(47958833331200, C[17]);
    try std.testing.expectEqual(10727718265600, C[18]);
    try std.testing.expectEqual(397478400000, C[19]);

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

    trmm(.row_major, .right, .lower, .no_trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(430, C[0]);
    try std.testing.expectEqual(456, C[1]);
    try std.testing.expectEqual(452, C[2]);
    try std.testing.expectEqual(392, C[3]);
    try std.testing.expectEqual(250, C[4]);
    try std.testing.expectEqual(980, C[5]);
    try std.testing.expectEqual(1036, C[6]);
    try std.testing.expectEqual(992, C[7]);
    try std.testing.expectEqual(822, C[8]);
    try std.testing.expectEqual(500, C[9]);
    try std.testing.expectEqual(1530, C[10]);
    try std.testing.expectEqual(1616, C[11]);
    try std.testing.expectEqual(1532, C[12]);
    try std.testing.expectEqual(1252, C[13]);
    try std.testing.expectEqual(750, C[14]);
    try std.testing.expectEqual(2080, C[15]);
    try std.testing.expectEqual(2196, C[16]);
    try std.testing.expectEqual(2072, C[17]);
    try std.testing.expectEqual(1682, C[18]);
    try std.testing.expectEqual(1000, C[19]);

    trmm(.col_major, .right, .lower, .no_trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(41008, C[0]);
    try std.testing.expectEqual(38568, C[1]);
    try std.testing.expectEqual(37048, C[2]);
    try std.testing.expectEqual(41088, C[3]);
    try std.testing.expectEqual(88148, C[4]);
    try std.testing.expectEqual(85696, C[5]);
    try std.testing.expectEqual(86124, C[6]);
    try std.testing.expectEqual(97184, C[7]);
    try std.testing.expectEqual(130148, C[8]);
    try std.testing.expectEqual(110216, C[9]);
    try std.testing.expectEqual(111240, C[10]);
    try std.testing.expectEqual(130256, C[11]);
    try std.testing.expectEqual(146056, C[12]);
    try std.testing.expectEqual(130456, C[13]);
    try std.testing.expectEqual(95780, C[14]);
    try std.testing.expectEqual(119040, C[15]);
    try std.testing.expectEqual(109800, C[16]);
    try std.testing.expectEqual(103600, C[17]);
    try std.testing.expectEqual(84100, C[18]);
    try std.testing.expectEqual(50000, C[19]);

    trmm(.row_major, .right, .lower, .no_trans, .unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(6376920, C[0]);
    try std.testing.expectEqual(6241792, C[1]);
    try std.testing.expectEqual(5608072, C[2]);
    try std.testing.expectEqual(4313280, C[3]);
    try std.testing.expectEqual(176296, C[4]);
    try std.testing.expectEqual(12136736, C[5]);
    try std.testing.expectEqual(11779200, C[6]);
    try std.testing.expectEqual(9949632, C[7]);
    try std.testing.expectEqual(5550664, C[8]);
    try std.testing.expectEqual(220432, C[9]);
    try std.testing.expectEqual(13196136, C[10]);
    try std.testing.expectEqual(12415680, C[11]);
    try std.testing.expectEqual(9394408, C[12]);
    try std.testing.expectEqual(4858352, C[13]);
    try std.testing.expectEqual(191560, C[14]);
    try std.testing.expectEqual(8626080, C[15]);
    try std.testing.expectEqual(7765400, C[16]);
    try std.testing.expectEqual(5534800, C[17]);
    try std.testing.expectEqual(2568200, C[18]);
    try std.testing.expectEqual(100000, C[19]);

    trmm(.col_major, .right, .lower, .no_trans, .unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(199572272, C[0]);
    try std.testing.expectEqual(156567936, C[1]);
    try std.testing.expectEqual(164724240, C[2]);
    try std.testing.expectEqual(192927808, C[3]);
    try std.testing.expectEqual(413570560, C[4]);
    try std.testing.expectEqual(225946720, C[5]);
    try std.testing.expectEqual(289508656, C[6]);
    try std.testing.expectEqual(375819584, C[7]);
    try std.testing.expectEqual(507106752, C[8]);
    try std.testing.expectEqual(302518720, C[9]);
    try std.testing.expectEqual(108801952, C[10]);
    try std.testing.expectEqual(269361600, C[11]);
    try std.testing.expectEqual(329404816, C[12]);
    try std.testing.expectEqual(231108704, C[13]);
    try std.testing.expectEqual(103111120, C[14]);
    try std.testing.expectEqual(21252160, C[15]);
    try std.testing.expectEqual(15530800, C[16]);
    try std.testing.expectEqual(11069600, C[17]);
    try std.testing.expectEqual(5136400, C[18]);
    try std.testing.expectEqual(200000, C[19]);

    trmm(.row_major, .right, .lower, .trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(399144544, C[0]);
    try std.testing.expectEqual(4586818368, C[1]);
    try std.testing.expectEqual(12431050688, C[2]);
    try std.testing.expectEqual(24970951872, C[3]);
    try std.testing.expectEqual(52787402432, C[4]);
    try std.testing.expectEqual(451893440, C[5]);
    try std.testing.expectEqual(6764481824, C[6]);
    try std.testing.expectEqual(21690344768, C[7]);
    try std.testing.expectEqual(49873150944, C[8]);
    try std.testing.expectEqual(78982904064, C[9]);
    try std.testing.expectEqual(217603904, C[10]);
    try std.testing.expectEqual(5076685824, C[11]);
    try std.testing.expectEqual(17422846560, C[12]);
    try std.testing.expectEqual(33280660992, C[13]);
    try std.testing.expectEqual(47822987712, C[14]);
    try std.testing.expectEqual(42504320, C[15]);
    try std.testing.expectEqual(472457120, C[16]);
    try std.testing.expectEqual(1128096320, C[17]);
    try std.testing.expectEqual(1801805120, C[18]);
    try std.testing.expectEqual(2341694720, C[19]);

    trmm(.col_major, .right, .lower, .trans, .non_unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(798289088, C[0]);
    try std.testing.expectEqual(9173636736, C[1]);
    try std.testing.expectEqual(24862101376, C[2]);
    try std.testing.expectEqual(49941903744, C[3]);
    try std.testing.expectEqual(740620212224, C[4]);
    try std.testing.expectEqual(24673781632, C[5]);
    try std.testing.expectEqual(144426948288, C[6]);
    try std.testing.expectEqual(403548634240, C[7]);
    try std.testing.expectEqual(2143695230720, C[8]);
    try std.testing.expectEqual(2088306710912, C[9]);
    try std.testing.expectEqual(188475714816, C[10]);
    try std.testing.expectEqual(628865058944, C[11]);
    try std.testing.expectEqual(3011882795840, C[12]);
    try std.testing.expectEqual(3521015060352, C[13]);
    try std.testing.expectEqual(2044575520704, C[14]);
    try std.testing.expectEqual(733956188032, C[15]);
    try std.testing.expectEqual(3276470740800, C[16]);
    try std.testing.expectEqual(3812024430080, C[17]);
    try std.testing.expectEqual(2269138024960, C[18]);
    try std.testing.expectEqual(954601897600, C[19]);

    trmm(.row_major, .right, .lower, .trans, .unit, m, n, alpha, B.ptr, n, C.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(1596578176, C[0]);
    try std.testing.expectEqual(27926742528, C[1]);
    try std.testing.expectEqual(287453844352, C[2]);
    try std.testing.expectEqual(1332368356864, C[3]);
    try std.testing.expectEqual(5459276625536, C[4]);
    try std.testing.expectEqual(49347563264, C[5]);
    try std.testing.expectEqual(584939276160, C[6]);
    try std.testing.expectEqual(4816167223296, C[7]);
    try std.testing.expectEqual(24515218548096, C[8]);
    try std.testing.expectEqual(133028306224640, C[9]);
    try std.testing.expectEqual(376951429632, C[10]);
    try std.testing.expectEqual(3519438695680, C[11]);
    try std.testing.expectEqual(25262992732288, C[12]);
    try std.testing.expectEqual(142882445649152, C[13]);
    try std.testing.expectEqual(347230525162752, C[14]);
    try std.testing.expectEqual(1467912376064, C[15]);
    try std.testing.expectEqual(15360415737984, C[16]);
    try std.testing.expectEqual(102406382776064, C[17]);
    try std.testing.expectEqual(276657758737024, C[18]);
    try std.testing.expectEqual(461171825269504, C[19]);

    trmm(.col_major, .right, .lower, .trans, .unit, m, n, alpha, B.ptr, n, C.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(3193156352, C[0]);
    try std.testing.expectEqual(55853485056, C[1]);
    try std.testing.expectEqual(574907688704, C[2]);
    try std.testing.expectEqual(2664736713728, C[3]);
    try std.testing.expectEqual(10924939563776, C[4]);
    try std.testing.expectEqual(210402096640, C[5]);
    try std.testing.expectEqual(2319693929728, C[6]);
    try std.testing.expectEqual(14961807874048, C[7]);
    try std.testing.expectEqual(136388442573824, C[8]);
    try std.testing.expectEqual(267013733916672, C[9]);
    try std.testing.expectEqual(11837654343936, C[10]);
    try std.testing.expectEqual(92091763105280, C[11]);
    try std.testing.expectEqual(835231856696320, C[12]);
    try std.testing.expectEqual(4011669135667200, C[13]);
    try std.testing.expectEqual(717844228080896, C[14]);
    try std.testing.expectEqual(198830065105408, C[15]);
    try std.testing.expectEqual(1885898595502848, C[16]);
    try std.testing.expectEqual(9912225996947968, C[17]);
    try std.testing.expectEqual(14468418390839808, C[18]);
    try std.testing.expectEqual(1196290334486528, C[19]);

    const beta = cf64.init(2, 2);

    const D = try a.alloc(cf64, m * m);
    defer a.free(D);
    const E = try a.alloc(cf64, n * n);
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
        cf64.init(13, 13),
        cf64.init(14, 14),
        cf64.init(15, 15),
        cf64.init(16, 16),
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

    trmm(.row_major, .left, .upper, .no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-440, F[0].re);
    try std.testing.expectEqual(440, F[0].im);
    try std.testing.expectEqual(-480, F[1].re);
    try std.testing.expectEqual(480, F[1].im);
    try std.testing.expectEqual(-520, F[2].re);
    try std.testing.expectEqual(520, F[2].im);
    try std.testing.expectEqual(-560, F[3].re);
    try std.testing.expectEqual(560, F[3].im);
    try std.testing.expectEqual(-600, F[4].re);
    try std.testing.expectEqual(600, F[4].im);
    try std.testing.expectEqual(-964, F[5].re);
    try std.testing.expectEqual(964, F[5].im);
    try std.testing.expectEqual(-1048, F[6].re);
    try std.testing.expectEqual(1048, F[6].im);
    try std.testing.expectEqual(-1132, F[7].re);
    try std.testing.expectEqual(1132, F[7].im);
    try std.testing.expectEqual(-1216, F[8].re);
    try std.testing.expectEqual(1216, F[8].im);
    try std.testing.expectEqual(-1300, F[9].re);
    try std.testing.expectEqual(1300, F[9].im);
    try std.testing.expectEqual(-1252, F[10].re);
    try std.testing.expectEqual(1252, F[10].im);
    try std.testing.expectEqual(-1344, F[11].re);
    try std.testing.expectEqual(1344, F[11].im);
    try std.testing.expectEqual(-1436, F[12].re);
    try std.testing.expectEqual(1436, F[12].im);
    try std.testing.expectEqual(-1528, F[13].re);
    try std.testing.expectEqual(1528, F[13].im);
    try std.testing.expectEqual(-1620, F[14].re);
    try std.testing.expectEqual(1620, F[14].im);
    try std.testing.expectEqual(-1024, F[15].re);
    try std.testing.expectEqual(1024, F[15].im);
    try std.testing.expectEqual(-1088, F[16].re);
    try std.testing.expectEqual(1088, F[16].im);
    try std.testing.expectEqual(-1152, F[17].re);
    try std.testing.expectEqual(1152, F[17].im);
    try std.testing.expectEqual(-1216, F[18].re);
    try std.testing.expectEqual(1216, F[18].im);
    try std.testing.expectEqual(-1280, F[19].re);
    try std.testing.expectEqual(1280, F[19].im);

    trmm(.col_major, .left, .upper, .no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-59200, F[0].re);
    try std.testing.expectEqual(-59200, F[0].im);
    try std.testing.expectEqual(-63680, F[1].re);
    try std.testing.expectEqual(-63680, F[1].im);
    try std.testing.expectEqual(-56480, F[2].re);
    try std.testing.expectEqual(-56480, F[2].im);
    try std.testing.expectEqual(-35840, F[3].re);
    try std.testing.expectEqual(-35840, F[3].im);
    try std.testing.expectEqual(-118272, F[4].re);
    try std.testing.expectEqual(-118272, F[4].im);
    try std.testing.expectEqual(-128448, F[5].re);
    try std.testing.expectEqual(-128448, F[5].im);
    try std.testing.expectEqual(-114032, F[6].re);
    try std.testing.expectEqual(-114032, F[6].im);
    try std.testing.expectEqual(-72448, F[7].re);
    try std.testing.expectEqual(-72448, F[7].im);
    try std.testing.expectEqual(-145824, F[8].re);
    try std.testing.expectEqual(-145824, F[8].im);
    try std.testing.expectEqual(-156544, F[9].re);
    try std.testing.expectEqual(-156544, F[9].im);
    try std.testing.expectEqual(-135728, F[10].re);
    try std.testing.expectEqual(-135728, F[10].im);
    try std.testing.expectEqual(-86016, F[11].re);
    try std.testing.expectEqual(-86016, F[11].im);
    try std.testing.expectEqual(-147872, F[12].re);
    try std.testing.expectEqual(-147872, F[12].im);
    try std.testing.expectEqual(-158816, F[13].re);
    try std.testing.expectEqual(-158816, F[13].im);
    try std.testing.expectEqual(-132720, F[14].re);
    try std.testing.expectEqual(-132720, F[14].im);
    try std.testing.expectEqual(-65536, F[15].re);
    try std.testing.expectEqual(-65536, F[15].im);
    try std.testing.expectEqual(-137728, F[16].re);
    try std.testing.expectEqual(-137728, F[16].im);
    try std.testing.expectEqual(-147968, F[17].re);
    try std.testing.expectEqual(-147968, F[17].im);
    try std.testing.expectEqual(-130304, F[18].re);
    try std.testing.expectEqual(-130304, F[18].im);
    try std.testing.expectEqual(-81920, F[19].re);
    try std.testing.expectEqual(-81920, F[19].im);

    trmm(.row_major, .left, .upper, .no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(3704896, F[0].re);
    try std.testing.expectEqual(-3941696, F[0].im);
    try std.testing.expectEqual(4148096, F[1].re);
    try std.testing.expectEqual(-4402816, F[1].im);
    try std.testing.expectEqual(4721536, F[2].re);
    try std.testing.expectEqual(-4947456, F[2].im);
    try std.testing.expectEqual(5157248, F[3].re);
    try std.testing.expectEqual(-5300608, F[3].im);
    try std.testing.expectEqual(4155712, F[4].re);
    try std.testing.expectEqual(-4628800, F[4].im);
    try std.testing.expectEqual(5897536, F[5].re);
    try std.testing.expectEqual(-6411328, F[5].im);
    try std.testing.expectEqual(6815744, F[6].re);
    try std.testing.expectEqual(-7271872, F[6].im);
    try std.testing.expectEqual(8875392, F[7].re);
    try std.testing.expectEqual(-9165184, F[7].im);
    try std.testing.expectEqual(8616576, F[8].re);
    try std.testing.expectEqual(-9199872, F[8].im);
    try std.testing.expectEqual(6337600, F[9].re);
    try std.testing.expectEqual(-6963776, F[9].im);
    try std.testing.expectEqual(3145728, F[10].re);
    try std.testing.expectEqual(-3688640, F[10].im);
    try std.testing.expectEqual(6610944, F[11].re);
    try std.testing.expectEqual(-6955008, F[11].im);
    try std.testing.expectEqual(7102464, F[12].re);
    try std.testing.expectEqual(-7693952, F[12].im);
    try std.testing.expectEqual(6254592, F[13].re);
    try std.testing.expectEqual(-6889856, F[13].im);
    try std.testing.expectEqual(3932160, F[14].re);
    try std.testing.expectEqual(-4463040, F[14].im);
    try std.testing.expectEqual(0, F[15].re);
    try std.testing.expectEqual(-262144, F[15].im);
    try std.testing.expectEqual(0, F[16].re);
    try std.testing.expectEqual(-550912, F[16].im);
    try std.testing.expectEqual(0, F[17].re);
    try std.testing.expectEqual(-591872, F[17].im);
    try std.testing.expectEqual(0, F[18].re);
    try std.testing.expectEqual(-521216, F[18].im);
    try std.testing.expectEqual(0, F[19].re);
    try std.testing.expectEqual(-327680, F[19].im);

    trmm(.col_major, .left, .upper, .no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(557089536, F[0].re);
    try std.testing.expectEqual(520640512, F[0].im);
    try std.testing.expectEqual(511834112, F[1].re);
    try std.testing.expectEqual(477157888, F[1].im);
    try std.testing.expectEqual(337374464, F[2].re);
    try std.testing.expectEqual(308983040, F[2].im);
    try std.testing.expectEqual(20915712, F[3].re);
    try std.testing.expectEqual(-286720, F[3].im);
    try std.testing.expectEqual(884172544, F[4].re);
    try std.testing.expectEqual(823891712, F[4].im);
    try std.testing.expectEqual(828742912, F[5].re);
    try std.testing.expectEqual(768624128, F[5].im);
    try std.testing.expectEqual(578086272, F[6].re);
    try std.testing.expectEqual(531611264, F[6].im);
    try std.testing.expectEqual(36081152, F[7].re);
    try std.testing.expectEqual(-579584, F[7].im);
    try std.testing.expectEqual(669359872, F[8].re);
    try std.testing.expectEqual(582600704, F[8].im);
    try std.testing.expectEqual(563628800, F[9].re);
    try std.testing.expectEqual(494789632, F[9].im);
    try std.testing.expectEqual(430969216, F[10].re);
    try std.testing.expectEqual(395570816, F[10].im);
    try std.testing.expectEqual(27131904, F[11].re);
    try std.testing.expectEqual(-688128, F[11].im);
    try std.testing.expectEqual(341690880, F[12].re);
    try std.testing.expectEqual(265466624, F[12].im);
    try std.testing.expectEqual(219490560, F[13].re);
    try std.testing.expectEqual(156015872, F[13].im);
    try std.testing.expectEqual(32519040, F[14].re);
    try std.testing.expectEqual(-1061760, F[14].im);
    try std.testing.expectEqual(524288, F[15].re);
    try std.testing.expectEqual(-524288, F[15].im);
    try std.testing.expectEqual(48742400, F[16].re);
    try std.testing.expectEqual(-1101824, F[16].im);
    try std.testing.expectEqual(40382464, F[17].re);
    try std.testing.expectEqual(-1183744, F[17].im);
    try std.testing.expectEqual(20703232, F[18].re);
    try std.testing.expectEqual(-1042432, F[18].im);
    try std.testing.expectEqual(655360, F[19].re);
    try std.testing.expectEqual(-655360, F[19].im);

    trmm(.row_major, .left, .upper, .conj_no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(14038320640, F[0].re);
    try std.testing.expectEqual(12970016256, F[0].im);
    try std.testing.expectEqual(7777487872, F[1].re);
    try std.testing.expectEqual(6135634944, F[1].im);
    try std.testing.expectEqual(6384557056, F[2].re);
    try std.testing.expectEqual(4397955072, F[2].im);
    try std.testing.expectEqual(8403680256, F[3].re);
    try std.testing.expectEqual(6515170304, F[3].im);
    try std.testing.expectEqual(8446434816, F[4].re);
    try std.testing.expectEqual(7230657024, F[4].im);
    try std.testing.expectEqual(31973745152, F[5].re);
    try std.testing.expectEqual(29506184704, F[5].im);
    try std.testing.expectEqual(16193520640, F[6].re);
    try std.testing.expectEqual(12704144384, F[6].im);
    try std.testing.expectEqual(11725531136, F[7].re);
    try std.testing.expectEqual(7381275648, F[7].im);
    try std.testing.expectEqual(22872876032, F[8].re);
    try std.testing.expectEqual(18317503488, F[8].im);
    try std.testing.expectEqual(14458595840, F[9].re);
    try std.testing.expectEqual(11824250368, F[9].im);
    try std.testing.expectEqual(18987811328, F[10].re);
    try std.testing.expectEqual(17379950080, F[10].im);
    try std.testing.expectEqual(3533438976, F[11].re);
    try std.testing.expectEqual(-83165184, F[11].im);
    try std.testing.expectEqual(16972756992, F[12].re);
    try std.testing.expectEqual(11623711744, F[12].im);
    try std.testing.expectEqual(10651339776, F[13].re);
    try std.testing.expectEqual(6814661632, F[13].im);
    try std.testing.expectEqual(1462295040, F[14].re);
    try std.testing.expectEqual(-78174720, F[14].im);
    try std.testing.expectEqual(33554432, F[15].re);
    try std.testing.expectEqual(-33554432, F[15].im);
    try std.testing.expectEqual(3119513600, F[16].re);
    try std.testing.expectEqual(-70516736, F[16].im);
    try std.testing.expectEqual(2584477696, F[17].re);
    try std.testing.expectEqual(-75759616, F[17].im);
    try std.testing.expectEqual(1325006848, F[18].re);
    try std.testing.expectEqual(-66715648, F[18].im);
    try std.testing.expectEqual(41943040, F[19].re);
    try std.testing.expectEqual(-41943040, F[19].im);

    trmm(.col_major, .left, .upper, .conj_no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(878538467328, F[0].re);
    try std.testing.expectEqual(671708002304, F[0].im);
    try std.testing.expectEqual(912648085504, F[1].re);
    try std.testing.expectEqual(688022978560, F[1].im);
    try std.testing.expectEqual(785141325824, F[2].re);
    try std.testing.expectEqual(584420241408, F[2].im);
    try std.testing.expectEqual(537835536384, F[3].re);
    try std.testing.expectEqual(416970899456, F[3].im);
    try std.testing.expectEqual(1865955004416, F[4].re);
    try std.testing.expectEqual(1460221853696, F[4].im);
    try std.testing.expectEqual(2071740452864, F[5].re);
    try std.testing.expectEqual(1629665644544, F[5].im);
    try std.testing.expectEqual(1416046776320, F[6].re);
    try std.testing.expectEqual(1001858891776, F[6].im);
    try std.testing.expectEqual(750433992704, F[7].re);
    try std.testing.expectEqual(472401641472, F[7].im);
    try std.testing.expectEqual(1247963455488, F[8].re);
    try std.testing.expectEqual(931108634624, F[8].im);
    try std.testing.expectEqual(1304391335936, F[9].re);
    try std.testing.expectEqual(974322761728, F[9].im);
    try std.testing.expectEqual(1047470036992, F[10].re);
    try std.testing.expectEqual(759727892480, F[10].im);
    try std.testing.expectEqual(226140094464, F[11].re);
    try std.testing.expectEqual(-5322571776, F[11].im);
    try std.testing.expectEqual(335305275392, F[12].re);
    try std.testing.expectEqual(178228959232, F[12].im);
    try std.testing.expectEqual(316003004416, F[13].re);
    try std.testing.expectEqual(158545842176, F[13].im);
    try std.testing.expectEqual(66354247680, F[14].re);
    try std.testing.expectEqual(-5452953600, F[14].im);
    try std.testing.expectEqual(2147483648, F[15].re);
    try std.testing.expectEqual(-2147483648, F[15].im);
    try std.testing.expectEqual(114048892928, F[16].re);
    try std.testing.expectEqual(-6380060672, F[16].im);
    try std.testing.expectEqual(117376548864, F[17].re);
    try std.testing.expectEqual(-6835666944, F[17].im);
    try std.testing.expectEqual(60816883712, F[18].re);
    try std.testing.expectEqual(-5452070912, F[18].im);
    try std.testing.expectEqual(2684354560, F[19].re);
    try std.testing.expectEqual(-2684354560, F[19].im);

    trmm(.row_major, .left, .upper, .conj_no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(29591584735232, F[0].re);
    try std.testing.expectEqual(25220193067008, F[0].im);
    try std.testing.expectEqual(16316087844864, F[1].re);
    try std.testing.expectEqual(11050261430272, F[1].im);
    try std.testing.expectEqual(12306602196992, F[2].re);
    try std.testing.expectEqual(8547713105920, F[2].im);
    try std.testing.expectEqual(14990543110144, F[3].re);
    try std.testing.expectEqual(11173798920192, F[3].im);
    try std.testing.expectEqual(12085797634048, F[4].re);
    try std.testing.expectEqual(14338550693888, F[4].im);
    try std.testing.expectEqual(30282030129152, F[5].re);
    try std.testing.expectEqual(28606473707520, F[5].im);
    try std.testing.expectEqual(10809862987776, F[6].re);
    try std.testing.expectEqual(4482617384960, F[6].im);
    try std.testing.expectEqual(13700661977088, F[7].re);
    try std.testing.expectEqual(7217340784640, F[7].im);
    try std.testing.expectEqual(11427934044160, F[8].re);
    try std.testing.expectEqual(8622961491968, F[8].im);
    try std.testing.expectEqual(2603955429376, F[9].re);
    try std.testing.expectEqual(4318846148608, F[9].im);
    try std.testing.expectEqual(678563504128, F[10].re);
    try std.testing.expectEqual(3511316643840, F[10].im);
    try std.testing.expectEqual(5937272193024, F[11].re);
    try std.testing.expectEqual(135392133120, F[11].im);
    try std.testing.expectEqual(5948226977792, F[12].re);
    try std.testing.expectEqual(698956455936, F[12].im);
    try std.testing.expectEqual(3234124742656, F[13].re);
    try std.testing.expectEqual(687398289408, F[13].im);
    try std.testing.expectEqual(272463421440, F[14].re);
    try std.testing.expectEqual(-7046430720, F[14].im);
    try std.testing.expectEqual(8589934592, F[15].re);
    try std.testing.expectEqual(0, F[15].im);
    try std.testing.expectEqual(240857907200, F[16].re);
    try std.testing.expectEqual(215337664512, F[16].im);
    try std.testing.expectEqual(248424431616, F[17].re);
    try std.testing.expectEqual(221081763840, F[17].im);
    try std.testing.expectEqual(132537909248, F[18].re);
    try std.testing.expectEqual(110729625600, F[18].im);
    try std.testing.expectEqual(10737418240, F[19].re);
    try std.testing.expectEqual(0, F[19].im);

    trmm(.col_major, .left, .upper, .conj_no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(1557610461052928, F[0].re);
    try std.testing.expectEqual(1219383999873024, F[0].im);
    try std.testing.expectEqual(1342266154876928, F[1].re);
    try std.testing.expectEqual(1022373962317824, F[1].im);
    try std.testing.expectEqual(906950364790784, F[2].re);
    try std.testing.expectEqual(712136565817344, F[2].im);
    try std.testing.expectEqual(7633488379904, F[3].re);
    try std.testing.expectEqual(52328684060672, F[3].im);
    try std.testing.expectEqual(1702724586831872, F[4].re);
    try std.testing.expectEqual(1161654117466112, F[4].im);
    try std.testing.expectEqual(1202982703071232, F[5].re);
    try std.testing.expectEqual(701252787011584, F[5].im);
    try std.testing.expectEqual(834694209830912, F[6].re);
    try std.testing.expectEqual(463625407823872, F[6].im);
    try std.testing.expectEqual(12966642384896, F[7].re);
    try std.testing.expectEqual(41836005523456, F[7].im);
    try std.testing.expectEqual(390855493877760, F[8].re);
    try std.testing.expectEqual(259926504144896, F[8].im);
    try std.testing.expectEqual(356200001536000, F[9].re);
    try std.testing.expectEqual(161880228364288, F[9].im);
    try std.testing.expectEqual(350570825302016, F[10].re);
    try std.testing.expectEqual(16503288283136, F[10].im);
    try std.testing.expectEqual(11603760119808, F[11].re);
    try std.testing.expectEqual(12145328652288, F[11].im);
    try std.testing.expectEqual(85436395667456, F[12].re);
    try std.testing.expectEqual(26788661149696, F[12].im);
    try std.testing.expectEqual(16473026101248, F[13].re);
    try std.testing.expectEqual(7561188835328, F[13].im);
    try std.testing.expectEqual(1074415779840, F[14].re);
    try std.testing.expectEqual(530833981440, F[14].im);
    try std.testing.expectEqual(17179869184, F[15].re);
    try std.testing.expectEqual(17179869184, F[15].im);
    try std.testing.expectEqual(10349239599104, F[16].re);
    try std.testing.expectEqual(9320292941824, F[16].im);
    try std.testing.expectEqual(5957497126912, F[17].re);
    try std.testing.expectEqual(5368197414912, F[17].im);
    try std.testing.expectEqual(687861661696, F[18].re);
    try std.testing.expectEqual(486535069696, F[18].im);
    try std.testing.expectEqual(21474836480, F[19].re);
    try std.testing.expectEqual(21474836480, F[19].im);

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

    trmm(.row_major, .left, .upper, .trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-4, F[0].re);
    try std.testing.expectEqual(4, F[0].im);
    try std.testing.expectEqual(-8, F[1].re);
    try std.testing.expectEqual(8, F[1].im);
    try std.testing.expectEqual(-12, F[2].re);
    try std.testing.expectEqual(12, F[2].im);
    try std.testing.expectEqual(-16, F[3].re);
    try std.testing.expectEqual(16, F[3].im);
    try std.testing.expectEqual(-20, F[4].re);
    try std.testing.expectEqual(20, F[4].im);
    try std.testing.expectEqual(-152, F[5].re);
    try std.testing.expectEqual(152, F[5].im);
    try std.testing.expectEqual(-184, F[6].re);
    try std.testing.expectEqual(184, F[6].im);
    try std.testing.expectEqual(-216, F[7].re);
    try std.testing.expectEqual(216, F[7].im);
    try std.testing.expectEqual(-248, F[8].re);
    try std.testing.expectEqual(248, F[8].im);
    try std.testing.expectEqual(-280, F[9].re);
    try std.testing.expectEqual(280, F[9].im);
    try std.testing.expectEqual(-664, F[10].re);
    try std.testing.expectEqual(664, F[10].im);
    try std.testing.expectEqual(-748, F[11].re);
    try std.testing.expectEqual(748, F[11].im);
    try std.testing.expectEqual(-832, F[12].re);
    try std.testing.expectEqual(832, F[12].im);
    try std.testing.expectEqual(-916, F[13].re);
    try std.testing.expectEqual(916, F[13].im);
    try std.testing.expectEqual(-1000, F[14].re);
    try std.testing.expectEqual(1000, F[14].im);
    try std.testing.expectEqual(-1760, F[15].re);
    try std.testing.expectEqual(1760, F[15].im);
    try std.testing.expectEqual(-1920, F[16].re);
    try std.testing.expectEqual(1920, F[16].im);
    try std.testing.expectEqual(-2080, F[17].re);
    try std.testing.expectEqual(2080, F[17].im);
    try std.testing.expectEqual(-2240, F[18].re);
    try std.testing.expectEqual(2240, F[18].im);
    try std.testing.expectEqual(-2400, F[19].re);
    try std.testing.expectEqual(2400, F[19].im);

    trmm(.col_major, .left, .upper, .trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-16, F[0].re);
    try std.testing.expectEqual(-16, F[0].im);
    try std.testing.expectEqual(-272, F[1].re);
    try std.testing.expectEqual(-272, F[1].im);
    try std.testing.expectEqual(-992, F[2].re);
    try std.testing.expectEqual(-992, F[2].im);
    try std.testing.expectEqual(-2400, F[3].re);
    try std.testing.expectEqual(-2400, F[3].im);
    try std.testing.expectEqual(-80, F[4].re);
    try std.testing.expectEqual(-80, F[4].im);
    try std.testing.expectEqual(-4048, F[5].re);
    try std.testing.expectEqual(-4048, F[5].im);
    try std.testing.expectEqual(-14896, F[6].re);
    try std.testing.expectEqual(-14896, F[6].im);
    try std.testing.expectEqual(-34416, F[7].re);
    try std.testing.expectEqual(-34416, F[7].im);
    try std.testing.expectEqual(-992, F[8].re);
    try std.testing.expectEqual(-992, F[8].im);
    try std.testing.expectEqual(-11680, F[9].re);
    try std.testing.expectEqual(-11680, F[9].im);
    try std.testing.expectEqual(-49344, F[10].re);
    try std.testing.expectEqual(-49344, F[10].im);
    try std.testing.expectEqual(-116288, F[11].re);
    try std.testing.expectEqual(-116288, F[11].im);
    try std.testing.expectEqual(-3328, F[12].re);
    try std.testing.expectEqual(-3328, F[12].im);
    try std.testing.expectEqual(-38624, F[13].re);
    try std.testing.expectEqual(-38624, F[13].im);
    try std.testing.expectEqual(-110592, F[14].re);
    try std.testing.expectEqual(-110592, F[14].im);
    try std.testing.expectEqual(-267200, F[15].re);
    try std.testing.expectEqual(-267200, F[15].im);
    try std.testing.expectEqual(-7680, F[16].re);
    try std.testing.expectEqual(-7680, F[16].im);
    try std.testing.expectEqual(-88320, F[17].re);
    try std.testing.expectEqual(-88320, F[17].im);
    try std.testing.expectEqual(-250880, F[18].re);
    try std.testing.expectEqual(-250880, F[18].im);
    try std.testing.expectEqual(-504320, F[19].re);
    try std.testing.expectEqual(-504320, F[19].im);

    trmm(.row_major, .left, .upper, .trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(0, F[0].re);
    try std.testing.expectEqual(-64, F[0].im);
    try std.testing.expectEqual(0, F[1].re);
    try std.testing.expectEqual(-1088, F[1].im);
    try std.testing.expectEqual(0, F[2].re);
    try std.testing.expectEqual(-3968, F[2].im);
    try std.testing.expectEqual(0, F[3].re);
    try std.testing.expectEqual(-9600, F[3].im);
    try std.testing.expectEqual(0, F[4].re);
    try std.testing.expectEqual(-320, F[4].im);
    try std.testing.expectEqual(128, F[5].re);
    try std.testing.expectEqual(-16320, F[5].im);
    try std.testing.expectEqual(2176, F[6].re);
    try std.testing.expectEqual(-61760, F[6].im);
    try std.testing.expectEqual(7936, F[7].re);
    try std.testing.expectEqual(-145600, F[7].im);
    try std.testing.expectEqual(19200, F[8].re);
    try std.testing.expectEqual(-23168, F[8].im);
    try std.testing.expectEqual(640, F[9].re);
    try std.testing.expectEqual(-47360, F[9].im);
    try std.testing.expectEqual(113536, F[10].re);
    try std.testing.expectEqual(-310912, F[10].im);
    try std.testing.expectEqual(420352, F[11].re);
    try std.testing.expectEqual(-885504, F[11].im);
    try std.testing.expectEqual(975552, F[12].re);
    try std.testing.expectEqual(-988864, F[12].im);
    try std.testing.expectEqual(56576, F[13].re);
    try std.testing.expectEqual(-211072, F[13].im);
    try std.testing.expectEqual(328000, F[14].re);
    try std.testing.expectEqual(-770368, F[14].im);
    try std.testing.expectEqual(2498304, F[15].re);
    try std.testing.expectEqual(-3567104, F[15].im);
    try std.testing.expectEqual(6062848, F[16].re);
    try std.testing.expectEqual(-6093568, F[16].im);
    try std.testing.expectEqual(1276928, F[17].re);
    try std.testing.expectEqual(-1630208, F[17].im);
    try std.testing.expectEqual(1924096, F[18].re);
    try std.testing.expectEqual(-2927616, F[18].im);
    try std.testing.expectEqual(5683456, F[19].re);
    try std.testing.expectEqual(-7700736, F[19].im);

    trmm(.col_major, .left, .upper, .trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(128, F[0].re);
    try std.testing.expectEqual(-128, F[0].im);
    try std.testing.expectEqual(3456, F[1].re);
    try std.testing.expectEqual(-2176, F[1].im);
    try std.testing.expectEqual(53760, F[2].re);
    try std.testing.expectEqual(-7936, F[2].im);
    try std.testing.expectEqual(321536, F[3].re);
    try std.testing.expectEqual(-19200, F[3].im);
    try std.testing.expectEqual(640, F[4].re);
    try std.testing.expectEqual(-640, F[4].im);
    try std.testing.expectEqual(39296, F[5].re);
    try std.testing.expectEqual(-32384, F[5].im);
    try std.testing.expectEqual(792192, F[6].re);
    try std.testing.expectEqual(-114048, F[6].im);
    try std.testing.expectEqual(4943232, F[7].re);
    try std.testing.expectEqual(-137600, F[7].im);
    try std.testing.expectEqual(84736, F[8].re);
    try std.testing.expectEqual(-7936, F[8].im);
    try std.testing.expectEqual(559360, F[9].re);
    try std.testing.expectEqual(290560, F[9].im);
    try std.testing.expectEqual(3577344, F[10].re);
    try std.testing.expectEqual(322048, F[10].im);
    try std.testing.expectEqual(25123328, F[11].re);
    try std.testing.expectEqual(6916096, F[11].im);
    try std.testing.expectEqual(3928832, F[12].re);
    try std.testing.expectEqual(-26624, F[12].im);
    try std.testing.expectEqual(20312576, F[13].re);
    try std.testing.expectEqual(19202048, F[13].im);
    try std.testing.expectEqual(46238720, F[14].re);
    try std.testing.expectEqual(36498176, F[14].im);
    try std.testing.expectEqual(121593856, F[15].re);
    try std.testing.expectEqual(71439360, F[15].im);
    try std.testing.expectEqual(24312832, F[16].re);
    try std.testing.expectEqual(-61440, F[16].im);
    try std.testing.expectEqual(127685632, F[17].re);
    try std.testing.expectEqual(120550400, F[17].im);
    try std.testing.expectEqual(294280192, F[18].re);
    try std.testing.expectEqual(267332608, F[18].im);
    try std.testing.expectEqual(610582528, F[19].re);
    try std.testing.expectEqual(498187264, F[19].im);

    trmm(.row_major, .left, .upper, .conj_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(512, F[0].re);
    try std.testing.expectEqual(-512, F[0].im);
    try std.testing.expectEqual(13824, F[1].re);
    try std.testing.expectEqual(-8704, F[1].im);
    try std.testing.expectEqual(215040, F[2].re);
    try std.testing.expectEqual(-31744, F[2].im);
    try std.testing.expectEqual(1286144, F[3].re);
    try std.testing.expectEqual(-76800, F[3].im);
    try std.testing.expectEqual(2560, F[4].re);
    try std.testing.expectEqual(-2560, F[4].im);
    try std.testing.expectEqual(944128, F[5].re);
    try std.testing.expectEqual(-778240, F[5].im);
    try std.testing.expectEqual(19040256, F[6].re);
    try std.testing.expectEqual(-2754560, F[6].im);
    try std.testing.expectEqual(119067648, F[7].re);
    try std.testing.expectEqual(-3365888, F[7].im);
    try std.testing.expectEqual(4605952, F[8].re);
    try std.testing.expectEqual(-344064, F[8].im);
    try std.testing.expectEqual(13429760, F[9].re);
    try std.testing.expectEqual(6968320, F[9].im);
    try std.testing.expectEqual(158504960, F[10].re);
    try std.testing.expectEqual(13261824, F[10].im);
    try std.testing.expectEqual(1127649280, F[11].re);
    try std.testing.expectEqual(301088768, F[11].im);
    try std.testing.expectEqual(311924224, F[12].re);
    try std.testing.expectEqual(-5119488, F[12].im);
    try std.testing.expectEqual(899984384, F[13].re);
    try std.testing.expectEqual(844437504, F[13].im);
    try std.testing.expectEqual(2050173440, F[14].re);
    try std.testing.expectEqual(1614047744, F[14].im);
    try std.testing.expectEqual(7954978816, F[15].re);
    try std.testing.expectEqual(4586539008, F[15].im);
    try std.testing.expectEqual(2787346432, F[16].re);
    try std.testing.expectEqual(324356096, F[16].im);
    try std.testing.expectEqual(8519507968, F[17].re);
    try std.testing.expectEqual(7709417472, F[17].im);
    try std.testing.expectEqual(19816792064, F[18].re);
    try std.testing.expectEqual(18030424064, F[18].im);
    try std.testing.expectEqual(41314650112, F[19].re);
    try std.testing.expectEqual(33645185024, F[19].im);

    trmm(.col_major, .left, .upper, .conj_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(2048, F[0].re);
    try std.testing.expectEqual(-2048, F[0].im);
    try std.testing.expectEqual(342016, F[1].re);
    try std.testing.expectEqual(-219136, F[1].im);
    try std.testing.expectEqual(10033152, F[2].re);
    try std.testing.expectEqual(-1763328, F[2].im);
    try std.testing.expectEqual(96016384, F[3].re);
    try std.testing.expectEqual(-7333888, F[3].im);
    try std.testing.expectEqual(10240, F[4].re);
    try std.testing.expectEqual(-10240, F[4].im);
    try std.testing.expectEqual(22710272, F[5].re);
    try std.testing.expectEqual(-18728960, F[5].im);
    try std.testing.expectEqual(875628544, F[6].re);
    try std.testing.expectEqual(-152422400, F[6].im);
    try std.testing.expectEqual(8815749120, F[7].re);
    try std.testing.expectEqual(-424404992, F[7].im);
    try std.testing.expectEqual(18423808, F[8].re);
    try std.testing.expectEqual(-1376256, F[8].im);
    try std.testing.expectEqual(414433280, F[9].re);
    try std.testing.expectEqual(160358400, F[9].im);
    try std.testing.expectEqual(7677222912, F[10].re);
    try std.testing.expectEqual(849866752, F[10].im);
    try std.testing.expectEqual(82671427584, F[11].re);
    try std.testing.expectEqual(20437725184, F[11].im);
    try std.testing.expectEqual(1247696896, F[12].re);
    try std.testing.expectEqual(-20477952, F[12].im);
    try std.testing.expectEqual(27838109696, F[13].re);
    try std.testing.expectEqual(20164110336, F[13].im);
    try std.testing.expectEqual(137436278784, F[14].re);
    try std.testing.expectEqual(104611299328, F[14].im);
    try std.testing.expectEqual(698748235776, F[15].re);
    try std.testing.expectEqual(437403648000, F[15].im);
    try std.testing.expectEqual(11149385728, F[16].re);
    try std.testing.expectEqual(1297424384, F[16].im);
    try std.testing.expectEqual(260215119872, F[17].re);
    try std.testing.expectEqual(191513141248, F[17].im);
    try std.testing.expectEqual(1313063641088, F[18].re);
    try std.testing.expectEqual(1113392177152, F[18].im);
    try std.testing.expectEqual(4455179591680, F[19].re);
    try std.testing.expectEqual(3683711180800, F[19].im);

    trmm(.row_major, .left, .upper, .conj_trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(8192, F[0].re);
    try std.testing.expectEqual(0, F[0].im);
    try std.testing.expectEqual(1122304, F[1].re);
    try std.testing.expectEqual(245760, F[1].im);
    try std.testing.expectEqual(23592960, F[2].re);
    try std.testing.expectEqual(16539648, F[2].im);
    try std.testing.expectEqual(206700544, F[3].re);
    try std.testing.expectEqual(177364992, F[3].im);
    try std.testing.expectEqual(40960, F[4].re);
    try std.testing.expectEqual(0, F[4].im);
    try std.testing.expectEqual(82894848, F[5].re);
    try std.testing.expectEqual(7946240, F[5].im);
    try std.testing.expectEqual(2058838016, F[6].re);
    try std.testing.expectEqual(1444659200, F[6].im);
    try std.testing.expectEqual(18560573440, F[7].re);
    try std.testing.expectEqual(16768581632, F[7].im);
    try std.testing.expectEqual(807731200, F[8].re);
    try std.testing.expectEqual(-24576000, F[8].im);
    try std.testing.expectEqual(508231680, F[9].re);
    try std.testing.expectEqual(1149501440, F[9].im);
    try std.testing.expectEqual(14290624512, F[10].re);
    try std.testing.expectEqual(16529743872, F[10].im);
    try std.testing.expectEqual(148989108224, F[11].re);
    try std.testing.expectEqual(201947848704, F[11].im);
    try std.testing.expectEqual(249497722880, F[12].re);
    try std.testing.expectEqual(-9450061824, F[12].im);
    try std.testing.expectEqual(17016061952, F[13].re);
    try std.testing.expectEqual(95877898240, F[13].im);
    try std.testing.expectEqual(77254213632, F[14].re);
    try std.testing.expectEqual(488585068544, F[14].im);
    try std.testing.expectEqual(891922636800, F[15].re);
    try std.testing.expectEqual(2312498012160, F[15].im);
    try std.testing.expectEqual(4015958032384, F[16].re);
    try std.testing.expectEqual(1001023406080, F[16].im);
    try std.testing.expectEqual(479557910528, F[17].re);
    try std.testing.expectEqual(888864407552, F[17].im);
    try std.testing.expectEqual(1737698017280, F[18].re);
    try std.testing.expectEqual(5820627550208, F[18].im);
    try std.testing.expectEqual(8153140232192, F[19].re);
    try std.testing.expectEqual(21304255217664, F[19].im);

    trmm(.col_major, .left, .upper, .conj_trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(16384, F[0].re);
    try std.testing.expectEqual(16384, F[0].im);
    try std.testing.expectEqual(1916928, F[1].re);
    try std.testing.expectEqual(2736128, F[1].im);
    try std.testing.expectEqual(59293696, F[2].re);
    try std.testing.expectEqual(90095616, F[2].im);
    try std.testing.expectEqual(1537523712, F[3].re);
    try std.testing.expectEqual(1774272512, F[3].im);
    try std.testing.expectEqual(81920, F[4].re);
    try std.testing.expectEqual(81920, F[4].im);
    try std.testing.expectEqual(150716416, F[5].re);
    try std.testing.expectEqual(181682176, F[5].im);
    try std.testing.expectEqual(4545626112, F[6].re);
    try std.testing.expectEqual(7324844032, F[6].im);
    try std.testing.expectEqual(131758505984, F[7].re);
    try std.testing.expectEqual(157782851584, F[7].im);
    try std.testing.expectEqual(1664614400, F[8].re);
    try std.testing.expectEqual(1566310400, F[8].im);
    try std.testing.expectEqual(14872084480, F[9].re);
    try std.testing.expectEqual(2823946240, F[9].im);
    try std.testing.expectEqual(44929351680, F[10].re);
    try std.testing.expectEqual(106736058368, F[10].im);
    try std.testing.expectEqual(821982986240, F[11].re);
    try std.testing.expectEqual(1756752674816, F[11].im);
    try std.testing.expectEqual(517895569408, F[12].re);
    try std.testing.expectEqual(480095322112, F[12].im);
    try std.testing.expectEqual(4832230785024, F[13].re);
    try std.testing.expectEqual(36786683904, F[13].im);
    try std.testing.expectEqual(8839898791936, F[14].re);
    try std.testing.expectEqual(4626592268288, F[14].im);
    try std.testing.expectEqual(15720883126272, F[15].re);
    try std.testing.expectEqual(40601704497152, F[15].im);
    try std.testing.expectEqual(6029869252608, F[16].re);
    try std.testing.expectEqual(10033962876928, F[16].im);
    try std.testing.expectEqual(79500547653632, F[17].re);
    try std.testing.expectEqual(22757312757760, F[17].im);
    try std.testing.expectEqual(155590946521088, F[18].re);
    try std.testing.expectEqual(86708070055936, F[18].im);
    try std.testing.expectEqual(313644711739392, F[19].re);
    try std.testing.expectEqual(509982067851264, F[19].im);

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

    trmm(.row_major, .left, .lower, .no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-4, F[0].re);
    try std.testing.expectEqual(4, F[0].im);
    try std.testing.expectEqual(-8, F[1].re);
    try std.testing.expectEqual(8, F[1].im);
    try std.testing.expectEqual(-12, F[2].re);
    try std.testing.expectEqual(12, F[2].im);
    try std.testing.expectEqual(-16, F[3].re);
    try std.testing.expectEqual(16, F[3].im);
    try std.testing.expectEqual(-20, F[4].re);
    try std.testing.expectEqual(20, F[4].im);
    try std.testing.expectEqual(-164, F[5].re);
    try std.testing.expectEqual(164, F[5].im);
    try std.testing.expectEqual(-208, F[6].re);
    try std.testing.expectEqual(208, F[6].im);
    try std.testing.expectEqual(-252, F[7].re);
    try std.testing.expectEqual(252, F[7].im);
    try std.testing.expectEqual(-296, F[8].re);
    try std.testing.expectEqual(296, F[8].im);
    try std.testing.expectEqual(-340, F[9].re);
    try std.testing.expectEqual(340, F[9].im);
    try std.testing.expectEqual(-760, F[10].re);
    try std.testing.expectEqual(760, F[10].im);
    try std.testing.expectEqual(-880, F[11].re);
    try std.testing.expectEqual(880, F[11].im);
    try std.testing.expectEqual(-1000, F[12].re);
    try std.testing.expectEqual(1000, F[12].im);
    try std.testing.expectEqual(-1120, F[13].re);
    try std.testing.expectEqual(1120, F[13].im);
    try std.testing.expectEqual(-1240, F[14].re);
    try std.testing.expectEqual(1240, F[14].im);
    try std.testing.expectEqual(-2072, F[15].re);
    try std.testing.expectEqual(2072, F[15].im);
    try std.testing.expectEqual(-2304, F[16].re);
    try std.testing.expectEqual(2304, F[16].im);
    try std.testing.expectEqual(-2536, F[17].re);
    try std.testing.expectEqual(2536, F[17].im);
    try std.testing.expectEqual(-2768, F[18].re);
    try std.testing.expectEqual(2768, F[18].im);
    try std.testing.expectEqual(-3000, F[19].re);
    try std.testing.expectEqual(3000, F[19].im);

    trmm(.col_major, .left, .lower, .no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-16, F[0].re);
    try std.testing.expectEqual(-16, F[0].im);
    try std.testing.expectEqual(-224, F[1].re);
    try std.testing.expectEqual(-224, F[1].im);
    try std.testing.expectEqual(-800, F[2].re);
    try std.testing.expectEqual(-800, F[2].im);
    try std.testing.expectEqual(-1920, F[3].re);
    try std.testing.expectEqual(-1920, F[3].im);
    try std.testing.expectEqual(-80, F[4].re);
    try std.testing.expectEqual(-80, F[4].im);
    try std.testing.expectEqual(-4096, F[5].re);
    try std.testing.expectEqual(-4096, F[5].im);
    try std.testing.expectEqual(-13984, F[6].re);
    try std.testing.expectEqual(-13984, F[6].im);
    try std.testing.expectEqual(-31680, F[7].re);
    try std.testing.expectEqual(-31680, F[7].im);
    try std.testing.expectEqual(-1184, F[8].re);
    try std.testing.expectEqual(-1184, F[8].im);
    try std.testing.expectEqual(-10528, F[9].re);
    try std.testing.expectEqual(-10528, F[9].im);
    try std.testing.expectEqual(-46512, F[10].re);
    try std.testing.expectEqual(-46512, F[10].im);
    try std.testing.expectEqual(-108416, F[11].re);
    try std.testing.expectEqual(-108416, F[11].im);
    try std.testing.expectEqual(-4000, F[12].re);
    try std.testing.expectEqual(-4000, F[12].im);
    try std.testing.expectEqual(-34880, F[13].re);
    try std.testing.expectEqual(-34880, F[13].im);
    try std.testing.expectEqual(-97920, F[14].re);
    try std.testing.expectEqual(-97920, F[14].im);
    try std.testing.expectEqual(-243968, F[15].re);
    try std.testing.expectEqual(-243968, F[15].im);
    try std.testing.expectEqual(-9216, F[16].re);
    try std.testing.expectEqual(-9216, F[16].im);
    try std.testing.expectEqual(-79296, F[17].re);
    try std.testing.expectEqual(-79296, F[17].im);
    try std.testing.expectEqual(-220448, F[18].re);
    try std.testing.expectEqual(-220448, F[18].im);
    try std.testing.expectEqual(-442880, F[19].re);
    try std.testing.expectEqual(-442880, F[19].im);

    trmm(.row_major, .left, .lower, .no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(0, F[0].re);
    try std.testing.expectEqual(-64, F[0].im);
    try std.testing.expectEqual(0, F[1].re);
    try std.testing.expectEqual(-896, F[1].im);
    try std.testing.expectEqual(0, F[2].re);
    try std.testing.expectEqual(-3200, F[2].im);
    try std.testing.expectEqual(0, F[3].re);
    try std.testing.expectEqual(-7680, F[3].im);
    try std.testing.expectEqual(0, F[4].re);
    try std.testing.expectEqual(-320, F[4].im);
    try std.testing.expectEqual(320, F[5].re);
    try std.testing.expectEqual(-16704, F[5].im);
    try std.testing.expectEqual(4480, F[6].re);
    try std.testing.expectEqual(-60416, F[6].im);
    try std.testing.expectEqual(16000, F[7].re);
    try std.testing.expectEqual(-142720, F[7].im);
    try std.testing.expectEqual(38400, F[8].re);
    try std.testing.expectEqual(-43136, F[8].im);
    try std.testing.expectEqual(1600, F[9].re);
    try std.testing.expectEqual(-43712, F[9].im);
    try std.testing.expectEqual(164416, F[10].re);
    try std.testing.expectEqual(-350464, F[10].im);
    try std.testing.expectEqual(567424, F[11].re);
    try std.testing.expectEqual(-1001088, F[11].im);
    try std.testing.expectEqual(1296000, F[12].re);
    try std.testing.expectEqual(-1312000, F[12].im);
    try std.testing.expectEqual(116480, F[13].re);
    try std.testing.expectEqual(-256000, F[13].im);
    try std.testing.expectEqual(424000, F[14].re);
    try std.testing.expectEqual(-815680, F[14].im);
    try std.testing.expectEqual(3020928, F[15].re);
    try std.testing.expectEqual(-3996800, F[15].im);
    try std.testing.expectEqual(7299712, F[16].re);
    try std.testing.expectEqual(-7336576, F[16].im);
    try std.testing.expectEqual(2055680, F[17].re);
    try std.testing.expectEqual(-2372864, F[17].im);
    try std.testing.expectEqual(2258944, F[18].re);
    try std.testing.expectEqual(-3140736, F[18].im);
    try std.testing.expectEqual(6468928, F[19].re);
    try std.testing.expectEqual(-8240448, F[19].im);

    trmm(.col_major, .left, .lower, .no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(128, F[0].re);
    try std.testing.expectEqual(-128, F[0].im);
    try std.testing.expectEqual(2304, F[1].re);
    try std.testing.expectEqual(-1792, F[1].im);
    try std.testing.expectEqual(32256, F[2].re);
    try std.testing.expectEqual(-6400, F[2].im);
    try std.testing.expectEqual(198656, F[3].re);
    try std.testing.expectEqual(-15360, F[3].im);
    try std.testing.expectEqual(640, F[4].re);
    try std.testing.expectEqual(-640, F[4].im);
    try std.testing.expectEqual(36608, F[5].re);
    try std.testing.expectEqual(-32768, F[5].im);
    try std.testing.expectEqual(601344, F[6].re);
    try std.testing.expectEqual(-102912, F[6].im);
    try std.testing.expectEqual(3757056, F[7].re);
    try std.testing.expectEqual(-28160, F[7].im);
    try std.testing.expectEqual(163072, F[8].re);
    try std.testing.expectEqual(-9472, F[8].im);
    try std.testing.expectEqual(435712, F[9].re);
    try std.testing.expectEqual(222976, F[9].im);
    try std.testing.expectEqual(2771328, F[10].re);
    try std.testing.expectEqual(133504, F[10].im);
    try std.testing.expectEqual(22048256, F[11].re);
    try std.testing.expectEqual(7690240, F[11].im);
    try std.testing.expectEqual(5216000, F[12].re);
    try std.testing.expectEqual(-32000, F[12].im);
    try std.testing.expectEqual(11240960, F[13].re);
    try std.testing.expectEqual(10088960, F[13].im);
    try std.testing.expectEqual(25391360, F[14].re);
    try std.testing.expectEqual(18030080, F[14].im);
    try std.testing.expectEqual(82372096, F[15].re);
    try std.testing.expectEqual(42863616, F[15].im);
    try std.testing.expectEqual(29272576, F[16].re);
    try std.testing.expectEqual(-73728, F[16].im);
    try std.testing.expectEqual(67549696, F[17].re);
    try std.testing.expectEqual(57763328, F[17].im);
    try std.testing.expectEqual(165278464, F[18].re);
    try std.testing.expectEqual(143392000, F[18].im);
    try std.testing.expectEqual(373490944, F[19].re);
    try std.testing.expectEqual(287463424, F[19].im);

    trmm(.row_major, .left, .lower, .conj_no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(512, F[0].re);
    try std.testing.expectEqual(-512, F[0].im);
    try std.testing.expectEqual(9216, F[1].re);
    try std.testing.expectEqual(-7168, F[1].im);
    try std.testing.expectEqual(129024, F[2].re);
    try std.testing.expectEqual(-25600, F[2].im);
    try std.testing.expectEqual(794624, F[3].re);
    try std.testing.expectEqual(-61440, F[3].im);
    try std.testing.expectEqual(2560, F[4].re);
    try std.testing.expectEqual(-2560, F[4].im);
    try std.testing.expectEqual(881152, F[5].re);
    try std.testing.expectEqual(-788992, F[5].im);
    try std.testing.expectEqual(14478336, F[6].re);
    try std.testing.expectEqual(-2505728, F[6].im);
    try std.testing.expectEqual(90814464, F[7].re);
    try std.testing.expectEqual(-803840, F[7].im);
    try std.testing.expectEqual(7886848, F[8].re);
    try std.testing.expectEqual(-534528, F[8].im);
    try std.testing.expectEqual(10469888, F[9].re);
    try std.testing.expectEqual(5338624, F[9].im);
    try std.testing.expectEqual(123407360, F[10].re);
    try std.testing.expectEqual(4558848, F[10].im);
    try std.testing.expectEqual(994259968, F[11].re);
    try std.testing.expectEqual(334189568, F[11].im);
    try std.testing.expectEqual(380947456, F[12].re);
    try std.testing.expectEqual(-2764800, F[12].im);
    try std.testing.expectEqual(508276736, F[13].re);
    try std.testing.expectEqual(442982400, F[13].im);
    try std.testing.expectEqual(1134671360, F[14].re);
    try std.testing.expectEqual(802219520, F[14].im);
    try std.testing.expectEqual(5440150528, F[15].re);
    try std.testing.expectEqual(2749440000, F[15].im);
    try std.testing.expectEqual(3230135296, F[16].re);
    try std.testing.expectEqual(450839552, F[16].im);
    try std.testing.expectEqual(4848212992, F[17].re);
    try std.testing.expectEqual(3693023232, F[17].im);
    try std.testing.expectEqual(11271741440, F[18].re);
    try std.testing.expectEqual(9781096448, F[18].im);
    try std.testing.expectEqual(25451335168, F[19].re);
    try std.testing.expectEqual(19491917312, F[19].im);

    trmm(.col_major, .left, .lower, .conj_no_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(2048, F[0].re);
    try std.testing.expectEqual(-2048, F[0].im);
    try std.testing.expectEqual(225280, F[1].re);
    try std.testing.expectEqual(-176128, F[1].im);
    try std.testing.expectEqual(5941248, F[2].re);
    try std.testing.expectEqual(-1333248, F[2].im);
    try std.testing.expectEqual(57352192, F[3].re);
    try std.testing.expectEqual(-5398528, F[3].im);
    try std.testing.expectEqual(10240, F[4].re);
    try std.testing.expectEqual(-10240, F[4].im);
    try std.testing.expectEqual(21168128, F[5].re);
    try std.testing.expectEqual(-18956288, F[5].im);
    try std.testing.expectEqual(661749760, F[6].re);
    try std.testing.expectEqual(-132374528, F[6].im);
    try std.testing.expectEqual(6535323648, F[7].re);
    try std.testing.expectEqual(-197009408, F[7].im);
    try std.testing.expectEqual(31547392, F[8].re);
    try std.testing.expectEqual(-2138112, F[8].im);
    try std.testing.expectEqual(314372096, F[9].re);
    try std.testing.expectEqual(123850752, F[9].im);
    try std.testing.expectEqual(5817722880, F[10].re);
    try std.testing.expectEqual(343656448, F[10].im);
    try std.testing.expectEqual(70017417216, F[11].re);
    try std.testing.expectEqual(21769240576, F[11].im);
    try std.testing.expectEqual(1523789824, F[12].re);
    try std.testing.expectEqual(-11059200, F[12].im);
    try std.testing.expectEqual(15246221312, F[13].re);
    try std.testing.expectEqual(10609459200, F[13].im);
    try std.testing.expectEqual(68728657920, F[14].re);
    try std.testing.expectEqual(47667988480, F[14].im);
    try std.testing.expectEqual(424993873920, F[15].re);
    try std.testing.expectEqual(228601896960, F[15].im);
    try std.testing.expectEqual(12920541184, F[16].re);
    try std.testing.expectEqual(1803358208, F[16].im);
    try std.testing.expectEqual(142198194176, F[17].re);
    try std.testing.expectEqual(92239273984, F[17].im);
    try std.testing.expectEqual(670468210688, F[18].re);
    try std.testing.expectEqual(539182968832, F[18].im);
    try std.testing.expectEqual(2376754020352, F[19].re);
    try std.testing.expectEqual(1842365513728, F[19].im);

    trmm(.row_major, .left, .lower, .conj_no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(8192, F[0].re);
    try std.testing.expectEqual(0, F[0].im);
    try std.testing.expectEqual(802816, F[1].re);
    try std.testing.expectEqual(98304, F[1].im);
    try std.testing.expectEqual(14548992, F[2].re);
    try std.testing.expectEqual(9216000, F[2].im);
    try std.testing.expectEqual(125501440, F[3].re);
    try std.testing.expectEqual(103907328, F[3].im);
    try std.testing.expectEqual(40960, F[4].re);
    try std.testing.expectEqual(0, F[4].im);
    try std.testing.expectEqual(80289792, F[5].re);
    try std.testing.expectEqual(4382720, F[5].im);
    try std.testing.expectEqual(1592754176, F[6].re);
    try std.testing.expectEqual(1055227904, F[6].im);
    try std.testing.expectEqual(13583491072, F[7].re);
    try std.testing.expectEqual(12649963520, F[7].im);
    try std.testing.expectEqual(1214414848, F[8].re);
    try std.testing.expectEqual(-49152000, F[8].im);
    try std.testing.expectEqual(381247488, F[9].re);
    try std.testing.expectEqual(876240896, F[9].im);
    try std.testing.expectEqual(11794931712, F[10].re);
    try std.testing.expectEqual(11564433408, F[10].im);
    try std.testing.expectEqual(122974453760, F[11].re);
    try std.testing.expectEqual(178271993856, F[11].im);
    try std.testing.expectEqual(264696528896, F[12].re);
    try std.testing.expectEqual(-4902912000, F[12].im);
    try std.testing.expectEqual(12600098816, F[13].re);
    try std.testing.expectEqual(51431489536, F[13].im);
    try std.testing.expectEqual(54696591360, F[14].re);
    try std.testing.expectEqual(237746954240, F[14].im);
    try std.testing.expectEqual(743032848384, F[15].re);
    try std.testing.expectEqual(1326749270016, F[15].im);
    try std.testing.expectEqual(4260349100032, F[16].re);
    try std.testing.expectEqual(1328180101120, F[16].im);
    try std.testing.expectEqual(557632299008, F[17].re);
    try std.testing.expectEqual(457109528576, F[17].im);
    try std.testing.expectEqual(1182092730368, F[18].re);
    try std.testing.expectEqual(3055469453312, F[18].im);
    try std.testing.expectEqual(5210101858304, F[19].re);
    try std.testing.expectEqual(11305253486592, F[19].im);

    trmm(.col_major, .left, .lower, .conj_no_trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(16384, F[0].re);
    try std.testing.expectEqual(16384, F[0].im);
    try std.testing.expectEqual(1474560, F[1].re);
    try std.testing.expectEqual(1802240, F[1].im);
    try std.testing.expectEqual(33243136, F[2].re);
    try std.testing.expectEqual(50282496, F[2].im);
    try std.testing.expectEqual(767361024, F[3].re);
    try std.testing.expectEqual(904331264, F[3].im);
    try std.testing.expectEqual(81920, F[4].re);
    try std.testing.expectEqual(81920, F[4].im);
    try std.testing.expectEqual(152141824, F[5].re);
    try std.testing.expectEqual(169345024, F[5].im);
    try std.testing.expectEqual(3323658240, F[6].re);
    try std.testing.expectEqual(5418680320, F[6].im);
    try std.testing.expectEqual(80889184256, F[7].re);
    try std.testing.expectEqual(103258095616, F[7].im);
    try std.testing.expectEqual(2527133696, F[8].re);
    try std.testing.expectEqual(2330525696, F[8].im);
    try std.testing.expectEqual(8725331968, F[9].re);
    try std.testing.expectEqual(2121760768, F[9].im);
    try std.testing.expectEqual(25708904448, F[10].re);
    try std.testing.expectEqual(70663651328, F[10].im);
    try std.testing.expectEqual(487192199168, F[11].re);
    try std.testing.expectEqual(1184838975488, F[11].im);
    try std.testing.expectEqual(539198881792, F[12].re);
    try std.testing.expectEqual(519587233792, F[12].im);
    try std.testing.expectEqual(2039909449728, F[13].re);
    try std.testing.expectEqual(88839880704, F[13].im);
    try std.testing.expectEqual(3163060387840, F[14].re);
    try std.testing.expectEqual(1966133854208, F[14].im);
    try std.testing.expectEqual(6096351166464, F[15].re);
    try std.testing.expectEqual(17118779113472, F[15].im);
    try std.testing.expectEqual(5864337997824, F[16].re);
    try std.testing.expectEqual(11177058402304, F[16].im);
    try std.testing.expectEqual(34283838341120, F[17].re);
    try std.testing.expectEqual(12654924464128, F[17].im);
    try std.testing.expectEqual(62991140126720, F[18].re);
    try std.testing.expectEqual(37212352380928, F[18].im);
    try std.testing.expectEqual(130559966969856, F[19].re);
    try std.testing.expectEqual(215571630981120, F[19].im);

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

    trmm(.row_major, .left, .lower, .trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-1352, F[0].re);
    try std.testing.expectEqual(1352, F[0].im);
    try std.testing.expectEqual(-1464, F[1].re);
    try std.testing.expectEqual(1464, F[1].im);
    try std.testing.expectEqual(-1576, F[2].re);
    try std.testing.expectEqual(1576, F[2].im);
    try std.testing.expectEqual(-1688, F[3].re);
    try std.testing.expectEqual(1688, F[3].im);
    try std.testing.expectEqual(-1800, F[4].re);
    try std.testing.expectEqual(1800, F[4].im);
    try std.testing.expectEqual(-1480, F[5].re);
    try std.testing.expectEqual(1480, F[5].im);
    try std.testing.expectEqual(-1600, F[6].re);
    try std.testing.expectEqual(1600, F[6].im);
    try std.testing.expectEqual(-1720, F[7].re);
    try std.testing.expectEqual(1720, F[7].im);
    try std.testing.expectEqual(-1840, F[8].re);
    try std.testing.expectEqual(1840, F[8].im);
    try std.testing.expectEqual(-1960, F[9].re);
    try std.testing.expectEqual(1960, F[9].im);
    try std.testing.expectEqual(-1444, F[10].re);
    try std.testing.expectEqual(1444, F[10].im);
    try std.testing.expectEqual(-1548, F[11].re);
    try std.testing.expectEqual(1548, F[11].im);
    try std.testing.expectEqual(-1652, F[12].re);
    try std.testing.expectEqual(1652, F[12].im);
    try std.testing.expectEqual(-1756, F[13].re);
    try std.testing.expectEqual(1756, F[13].im);
    try std.testing.expectEqual(-1860, F[14].re);
    try std.testing.expectEqual(1860, F[14].im);
    try std.testing.expectEqual(-1024, F[15].re);
    try std.testing.expectEqual(1024, F[15].im);
    try std.testing.expectEqual(-1088, F[16].re);
    try std.testing.expectEqual(1088, F[16].im);
    try std.testing.expectEqual(-1152, F[17].re);
    try std.testing.expectEqual(1152, F[17].im);
    try std.testing.expectEqual(-1216, F[18].re);
    try std.testing.expectEqual(1216, F[18].im);
    try std.testing.expectEqual(-1280, F[19].re);
    try std.testing.expectEqual(1280, F[19].im);

    trmm(.col_major, .left, .lower, .trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-63040, F[0].re);
    try std.testing.expectEqual(-63040, F[0].im);
    try std.testing.expectEqual(-133280, F[1].re);
    try std.testing.expectEqual(-133280, F[1].im);
    try std.testing.expectEqual(-150368, F[2].re);
    try std.testing.expectEqual(-150368, F[2].im);
    try std.testing.expectEqual(-108032, F[3].re);
    try std.testing.expectEqual(-108032, F[3].im);
    try std.testing.expectEqual(-65760, F[4].re);
    try std.testing.expectEqual(-65760, F[4].im);
    try std.testing.expectEqual(-135360, F[5].re);
    try std.testing.expectEqual(-135360, F[5].im);
    try std.testing.expectEqual(-152960, F[6].re);
    try std.testing.expectEqual(-152960, F[6].im);
    try std.testing.expectEqual(-110080, F[7].re);
    try std.testing.expectEqual(-110080, F[7].im);
    try std.testing.expectEqual(-65136, F[8].re);
    try std.testing.expectEqual(-65136, F[8].im);
    try std.testing.expectEqual(-137008, F[9].re);
    try std.testing.expectEqual(-137008, F[9].im);
    try std.testing.expectEqual(-137840, F[10].re);
    try std.testing.expectEqual(-137840, F[10].im);
    try std.testing.expectEqual(-99072, F[11].re);
    try std.testing.expectEqual(-99072, F[11].im);
    try std.testing.expectEqual(-59360, F[12].re);
    try std.testing.expectEqual(-59360, F[12].im);
    try std.testing.expectEqual(-126992, F[13].re);
    try std.testing.expectEqual(-126992, F[13].im);
    try std.testing.expectEqual(-130992, F[14].re);
    try std.testing.expectEqual(-130992, F[14].im);
    try std.testing.expectEqual(-65536, F[15].re);
    try std.testing.expectEqual(-65536, F[15].im);
    try std.testing.expectEqual(-48640, F[16].re);
    try std.testing.expectEqual(-48640, F[16].im);
    try std.testing.expectEqual(-102656, F[17].re);
    try std.testing.expectEqual(-102656, F[17].im);
    try std.testing.expectEqual(-114944, F[18].re);
    try std.testing.expectEqual(-114944, F[18].im);
    try std.testing.expectEqual(-81920, F[19].re);
    try std.testing.expectEqual(-81920, F[19].im);

    trmm(.row_major, .left, .lower, .trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(11077312, F[0].re);
    try std.testing.expectEqual(-11329472, F[0].im);
    try std.testing.expectEqual(9155072, F[1].re);
    try std.testing.expectEqual(-9688192, F[1].im);
    try std.testing.expectEqual(9676672, F[2].re);
    try std.testing.expectEqual(-10278144, F[2].im);
    try std.testing.expectEqual(11851520, F[3].re);
    try std.testing.expectEqual(-12283648, F[3].im);
    try std.testing.expectEqual(11715712, F[4].re);
    try std.testing.expectEqual(-11978752, F[4].im);
    try std.testing.expectEqual(9183616, F[5].re);
    try std.testing.expectEqual(-9725056, F[5].im);
    try std.testing.expectEqual(6686720, F[6].re);
    try std.testing.expectEqual(-7298560, F[6].im);
    try std.testing.expectEqual(8123136, F[7].re);
    try std.testing.expectEqual(-8563456, F[7].im);
    try std.testing.expectEqual(11516544, F[8].re);
    try std.testing.expectEqual(-11777088, F[8].im);
    try std.testing.expectEqual(9827200, F[9].re);
    try std.testing.expectEqual(-10375232, F[9].im);
    try std.testing.expectEqual(3932160, F[10].re);
    try std.testing.expectEqual(-4483520, F[10].im);
    try std.testing.expectEqual(2918400, F[11].re);
    try std.testing.expectEqual(-3314688, F[11].im);
    try std.testing.expectEqual(6159360, F[12].re);
    try std.testing.expectEqual(-6396800, F[12].im);
    try std.testing.expectEqual(6896640, F[13].re);
    try std.testing.expectEqual(-7404608, F[13].im);
    try std.testing.expectEqual(4915200, F[14].re);
    try std.testing.expectEqual(-5439168, F[14].im);
    try std.testing.expectEqual(0, F[15].re);
    try std.testing.expectEqual(-262144, F[15].im);
    try std.testing.expectEqual(0, F[16].re);
    try std.testing.expectEqual(-194560, F[16].im);
    try std.testing.expectEqual(0, F[17].re);
    try std.testing.expectEqual(-410624, F[17].im);
    try std.testing.expectEqual(0, F[18].re);
    try std.testing.expectEqual(-459776, F[18].im);
    try std.testing.expectEqual(0, F[19].re);
    try std.testing.expectEqual(-327680, F[19].im);

    trmm(.col_major, .left, .lower, .trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(442195200, F[0].re);
    try std.testing.expectEqual(378480640, F[0].im);
    try std.testing.expectEqual(718551296, F[1].re);
    try std.testing.expectEqual(649129216, F[1].im);
    try std.testing.expectEqual(629524736, F[2].re);
    try std.testing.expectEqual(567670016, F[2].im);
    try std.testing.expectEqual(48270336, F[3].re);
    try std.testing.expectEqual(-864256, F[3].im);
    try std.testing.expectEqual(349787392, F[4].re);
    try std.testing.expectEqual(283153664, F[4].im);
    try std.testing.expectEqual(516207616, F[5].re);
    try std.testing.expectEqual(446085632, F[5].im);
    try std.testing.expectEqual(439016448, F[6].re);
    try std.testing.expectEqual(388686848, F[6].im);
    try std.testing.expectEqual(33373184, F[7].re);
    try std.testing.expectEqual(-880640, F[7].im);
    try std.testing.expectEqual(236426368, F[8].re);
    try std.testing.expectEqual(171976832, F[8].im);
    try std.testing.expectEqual(272013440, F[9].re);
    try std.testing.expectEqual(202393216, F[9].im);
    try std.testing.expectEqual(175936384, F[10].re);
    try std.testing.expectEqual(138980480, F[10].im);
    try std.testing.expectEqual(12466176, F[11].re);
    try std.testing.expectEqual(-792576, F[11].im);
    try std.testing.expectEqual(153813504, F[12].re);
    try std.testing.expectEqual(113680640, F[12].im);
    try std.testing.expectEqual(189287808, F[13].re);
    try std.testing.expectEqual(136609664, F[13].im);
    try std.testing.expectEqual(33291648, F[14].re);
    try std.testing.expectEqual(-1047936, F[14].im);
    try std.testing.expectEqual(524288, F[15].re);
    try std.testing.expectEqual(-524288, F[15].im);
    try std.testing.expectEqual(14434304, F[16].re);
    try std.testing.expectEqual(-389120, F[16].im);
    try std.testing.expectEqual(24180736, F[17].re);
    try std.testing.expectEqual(-821248, F[17].im);
    try std.testing.expectEqual(16648192, F[18].re);
    try std.testing.expectEqual(-919552, F[18].im);
    try std.testing.expectEqual(655360, F[19].re);
    try std.testing.expectEqual(-655360, F[19].im);

    trmm(.row_major, .left, .lower, .conj_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(18453905920, F[0].re);
    try std.testing.expectEqual(15411669504, F[0].im);
    try std.testing.expectEqual(12853900288, F[1].re);
    try std.testing.expectEqual(10321486848, F[1].im);
    try std.testing.expectEqual(9980247040, F[2].re);
    try std.testing.expectEqual(6302865408, F[2].im);
    try std.testing.expectEqual(12601675776, F[3].re);
    try std.testing.expectEqual(8306210816, F[3].im);
    try std.testing.expectEqual(8071996416, F[4].re);
    try std.testing.expectEqual(5108674560, F[4].im);
    try std.testing.expectEqual(19455798272, F[5].re);
    try std.testing.expectEqual(16235914240, F[5].im);
    try std.testing.expectEqual(11843362816, F[6].re);
    try std.testing.expectEqual(9274990592, F[6].im);
    try std.testing.expectEqual(8307617792, F[7].re);
    try std.testing.expectEqual(4480100352, F[7].im);
    try std.testing.expectEqual(14178043904, F[8].re);
    try std.testing.expectEqual(9540335616, F[8].im);
    try std.testing.expectEqual(7896688640, F[9].re);
    try std.testing.expectEqual(4778819584, F[9].im);
    try std.testing.expectEqual(7772658176, F[10].re);
    try std.testing.expectEqual(6083683840, F[10].im);
    try std.testing.expectEqual(1414569984, F[11].re);
    try std.testing.expectEqual(-58220544, F[11].im);
    try std.testing.expectEqual(8218638336, F[12].re);
    try std.testing.expectEqual(4952673280, F[12].im);
    try std.testing.expectEqual(9327555072, F[13].re);
    try std.testing.expectEqual(5955652096, F[13].im);
    try std.testing.expectEqual(1504154112, F[14].re);
    try std.testing.expectEqual(-85430784, F[14].im);
    try std.testing.expectEqual(33554432, F[15].re);
    try std.testing.expectEqual(-33554432, F[15].im);
    try std.testing.expectEqual(923795456, F[16].re);
    try std.testing.expectEqual(-24903680, F[16].im);
    try std.testing.expectEqual(1547567104, F[17].re);
    try std.testing.expectEqual(-52559872, F[17].im);
    try std.testing.expectEqual(1065484288, F[18].re);
    try std.testing.expectEqual(-58851328, F[18].im);
    try std.testing.expectEqual(41943040, F[19].re);
    try std.testing.expectEqual(-41943040, F[19].im);

    trmm(.col_major, .left, .lower, .conj_trans, .non_unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(498036602880, F[0].re);
    try std.testing.expectEqual(352752330752, F[0].im);
    try std.testing.expectEqual(991194148864, F[1].re);
    try std.testing.expectEqual(689994661888, F[1].im);
    try std.testing.expectEqual(1044011307008, F[2].re);
    try std.testing.expectEqual(676024197120, F[2].im);
    try std.testing.expectEqual(806507249664, F[3].re);
    try std.testing.expectEqual(531597492224, F[3].im);
    try std.testing.expectEqual(462976610304, F[4].re);
    try std.testing.expectEqual(333303504896, F[4].im);
    try std.testing.expectEqual(1064397086720, F[5].re);
    try std.testing.expectEqual(792724889600, F[5].im);
    try std.testing.expectEqual(919873617920, F[6].re);
    try std.testing.expectEqual(623144402944, F[6].im);
    try std.testing.expectEqual(531687538688, F[7].re);
    try std.testing.expectEqual(286726422528, F[7].im);
    try std.testing.expectEqual(235790702592, F[8].re);
    try std.testing.expectEqual(148464576512, F[8].im);
    try std.testing.expectEqual(452421195776, F[9].re);
    try std.testing.expectEqual(283171760128, F[9].im);
    try std.testing.expectEqual(409896318976, F[10].re);
    try std.testing.expectEqual(264887502848, F[10].im);
    try std.testing.expectEqual(90532478976, F[11].re);
    try std.testing.expectEqual(-3726114816, F[11].im);
    try std.testing.expectEqual(126081714176, F[12].re);
    try std.testing.expectEqual(65893869568, F[12].im);
    try std.testing.expectEqual(267051378688, F[13].re);
    try std.testing.expectEqual(139469846528, F[13].im);
    try std.testing.expectEqual(67793393664, F[14].re);
    try std.testing.expectEqual(-5369567232, F[14].im);
    try std.testing.expectEqual(2147483648, F[15].re);
    try std.testing.expectEqual(-2147483648, F[15].im);
    try std.testing.expectEqual(29532618752, F[16].re);
    try std.testing.expectEqual(-1897398272, F[16].im);
    try std.testing.expectEqual(68317347840, F[17].re);
    try std.testing.expectEqual(-4251451392, F[17].im);
    try std.testing.expectEqual(48894574592, F[18].re);
    try std.testing.expectEqual(-4602724352, F[18].im);
    try std.testing.expectEqual(2684354560, F[19].re);
    try std.testing.expectEqual(-2684354560, F[19].im);

    trmm(.row_major, .left, .lower, .conj_trans, .unit, m, n, beta, D.ptr, m, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(36446446911488, F[0].re);
    try std.testing.expectEqual(26980356612096, F[0].im);
    try std.testing.expectEqual(23794736750592, F[1].re);
    try std.testing.expectEqual(15592460836864, F[1].im);
    try std.testing.expectEqual(19461168791552, F[2].re);
    try std.testing.expectEqual(11325703290880, F[2].im);
    try std.testing.expectEqual(17422001078272, F[3].re);
    try std.testing.expectEqual(10427073822720, F[3].im);
    try std.testing.expectEqual(11887918735360, F[4].re);
    try std.testing.expectEqual(6923104575488, F[4].im);
    try std.testing.expectEqual(17059456237568, F[5].re);
    try std.testing.expectEqual(14189484982272, F[5].im);
    try std.testing.expectEqual(5868584239104, F[6].re);
    try std.testing.expectEqual(2830737145856, F[6].im);
    try std.testing.expectEqual(9358962278400, F[7].re);
    try std.testing.expectEqual(4034501427200, F[7].im);
    try std.testing.expectEqual(13594803576832, F[8].re);
    try std.testing.expectEqual(6089551855616, F[8].im);
    try std.testing.expectEqual(3200558473216, F[9].re);
    try std.testing.expectEqual(1106079367168, F[9].im);
    try std.testing.expectEqual(418866651136, F[10].re);
    try std.testing.expectEqual(1220718624768, F[10].im);
    try std.testing.expectEqual(1960474312704, F[11].re);
    try std.testing.expectEqual(59768832000, F[11].im);
    try std.testing.expectEqual(4219416559616, F[12].re);
    try std.testing.expectEqual(128864083968, F[12].im);
    try std.testing.expectEqual(3188837539840, F[13].re);
    try std.testing.expectEqual(536878989312, F[13].im);
    try std.testing.expectEqual(307387195392, F[14].re);
    try std.testing.expectEqual(-36213620736, F[14].im);
    try std.testing.expectEqual(8589934592, F[15].re);
    try std.testing.expectEqual(0, F[15].im);
    try std.testing.expectEqual(62860034048, F[16].re);
    try std.testing.expectEqual(55270440960, F[16].im);
    try std.testing.expectEqual(145137598464, F[17].re);
    try std.testing.expectEqual(128131792896, F[17].im);
    try std.testing.expectEqual(106994597888, F[18].re);
    try std.testing.expectEqual(88583700480, F[18].im);
    try std.testing.expectEqual(10737418240, F[19].re);
    try std.testing.expectEqual(0, F[19].im);

    trmm(.col_major, .left, .lower, .conj_trans, .unit, m, n, beta, D.ptr, m, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(721576117354496, F[0].re);
    try std.testing.expectEqual(554334914396160, F[0].im);
    try std.testing.expectEqual(1118821312495616, F[1].re);
    try std.testing.expectEqual(729560449646592, F[1].im);
    try std.testing.expectEqual(852526982758400, F[2].re);
    try std.testing.expectEqual(562073287655424, F[2].im);
    try std.testing.expectEqual(13989854511104, F[3].re);
    try std.testing.expectEqual(55698149801984, F[3].im);
    try std.testing.expectEqual(366571685543936, F[4].re);
    try std.testing.expectEqual(249658795065344, F[4].im);
    try std.testing.expectEqual(469547094114304, F[5].re);
    try std.testing.expectEqual(270862568194048, F[5].im);
    try std.testing.expectEqual(455305883549696, F[6].re);
    try std.testing.expectEqual(211054711275520, F[6].im);
    try std.testing.expectEqual(10648921702400, F[7].re);
    try std.testing.expectEqual(26786927411200, F[7].im);
    try std.testing.expectEqual(77008960045056, F[8].re);
    try std.testing.expectEqual(63822270611456, F[8].im);
    try std.testing.expectEqual(78652402450432, F[9].re);
    try std.testing.expectEqual(44705999798272, F[9].im);
    try std.testing.expectEqual(92499063062528, F[10].re);
    try std.testing.expectEqual(6148074487808, F[10].im);
    try std.testing.expectEqual(3801410961408, F[11].re);
    try std.testing.expectEqual(4040486289408, F[11].im);
    try std.testing.expectEqual(37517890568192, F[12].re);
    try std.testing.expectEqual(12557029752832, F[12].im);
    try std.testing.expectEqual(14185636478976, F[13].re);
    try std.testing.expectEqual(6437451677696, F[13].im);
    try std.testing.expectEqual(1099518492672, F[14].re);
    try std.testing.expectEqual(542347149312, F[14].im);
    try std.testing.expectEqual(17179869184, F[15].re);
    try std.testing.expectEqual(17179869184, F[15].im);
    try std.testing.expectEqual(2632013840384, F[16].re);
    try std.testing.expectEqual(2324319698944, F[16].im);
    try std.testing.expectEqual(3373457735680, F[17].re);
    try std.testing.expectEqual(3026882396160, F[17].im);
    try std.testing.expectEqual(552217870336, F[18].re);
    try std.testing.expectEqual(391156596736, F[18].im);
    try std.testing.expectEqual(21474836480, F[19].re);
    try std.testing.expectEqual(21474836480, F[19].im);

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

    trmm(.row_major, .right, .upper, .no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-4, F[0].re);
    try std.testing.expectEqual(4, F[0].im);
    try std.testing.expectEqual(-64, F[1].re);
    try std.testing.expectEqual(64, F[1].im);
    try std.testing.expectEqual(-232, F[2].re);
    try std.testing.expectEqual(232, F[2].im);
    try std.testing.expectEqual(-560, F[3].re);
    try std.testing.expectEqual(560, F[3].im);
    try std.testing.expectEqual(-1100, F[4].re);
    try std.testing.expectEqual(1100, F[4].im);
    try std.testing.expectEqual(-24, F[5].re);
    try std.testing.expectEqual(24, F[5].im);
    try std.testing.expectEqual(-244, F[6].re);
    try std.testing.expectEqual(244, F[6].im);
    try std.testing.expectEqual(-712, F[7].re);
    try std.testing.expectEqual(712, F[7].im);
    try std.testing.expectEqual(-1480, F[8].re);
    try std.testing.expectEqual(1480, F[8].im);
    try std.testing.expectEqual(-2600, F[9].re);
    try std.testing.expectEqual(2600, F[9].im);
    try std.testing.expectEqual(-44, F[10].re);
    try std.testing.expectEqual(44, F[10].im);
    try std.testing.expectEqual(-424, F[11].re);
    try std.testing.expectEqual(424, F[11].im);
    try std.testing.expectEqual(-1192, F[12].re);
    try std.testing.expectEqual(1192, F[12].im);
    try std.testing.expectEqual(-2400, F[13].re);
    try std.testing.expectEqual(2400, F[13].im);
    try std.testing.expectEqual(-4100, F[14].re);
    try std.testing.expectEqual(4100, F[14].im);
    try std.testing.expectEqual(-64, F[15].re);
    try std.testing.expectEqual(64, F[15].im);
    try std.testing.expectEqual(-604, F[16].re);
    try std.testing.expectEqual(604, F[16].im);
    try std.testing.expectEqual(-1672, F[17].re);
    try std.testing.expectEqual(1672, F[17].im);
    try std.testing.expectEqual(-3320, F[18].re);
    try std.testing.expectEqual(3320, F[18].im);
    try std.testing.expectEqual(-5600, F[19].re);
    try std.testing.expectEqual(5600, F[19].im);

    trmm(.col_major, .right, .upper, .no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-16, F[0].re);
    try std.testing.expectEqual(-16, F[0].im);
    try std.testing.expectEqual(-256, F[1].re);
    try std.testing.expectEqual(-256, F[1].im);
    try std.testing.expectEqual(-928, F[2].re);
    try std.testing.expectEqual(-928, F[2].im);
    try std.testing.expectEqual(-2240, F[3].re);
    try std.testing.expectEqual(-2240, F[3].im);
    try std.testing.expectEqual(-30896, F[4].re);
    try std.testing.expectEqual(-30896, F[4].im);
    try std.testing.expectEqual(-2208, F[5].re);
    try std.testing.expectEqual(-2208, F[5].im);
    try std.testing.expectEqual(-12400, F[6].re);
    try std.testing.expectEqual(-12400, F[6].im);
    try std.testing.expectEqual(-33376, F[7].re);
    try std.testing.expectEqual(-33376, F[7].im);
    try std.testing.expectEqual(-129936, F[8].re);
    try std.testing.expectEqual(-129936, F[8].im);
    try std.testing.expectEqual(-139168, F[9].re);
    try std.testing.expectEqual(-139168, F[9].im);
    try std.testing.expectEqual(-24208, F[10].re);
    try std.testing.expectEqual(-24208, F[10].im);
    try std.testing.expectEqual(-80864, F[11].re);
    try std.testing.expectEqual(-80864, F[11].im);
    try std.testing.expectEqual(-272208, F[12].re);
    try std.testing.expectEqual(-272208, F[12].im);
    try std.testing.expectEqual(-375328, F[13].re);
    try std.testing.expectEqual(-375328, F[13].im);
    try std.testing.expectEqual(-346208, F[14].re);
    try std.testing.expectEqual(-346208, F[14].im);
    try std.testing.expectEqual(-119648, F[15].re);
    try std.testing.expectEqual(-119648, F[15].im);
    try std.testing.expectEqual(-408128, F[16].re);
    try std.testing.expectEqual(-408128, F[16].im);
    try std.testing.expectEqual(-644288, F[17].re);
    try std.testing.expectEqual(-644288, F[17].im);
    try std.testing.expectEqual(-770608, F[18].re);
    try std.testing.expectEqual(-770608, F[18].im);
    try std.testing.expectEqual(-714848, F[19].re);
    try std.testing.expectEqual(-714848, F[19].im);

    trmm(.row_major, .right, .upper, .no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(0, F[0].re);
    try std.testing.expectEqual(-64, F[0].im);
    try std.testing.expectEqual(128, F[1].re);
    try std.testing.expectEqual(-1152, F[1].im);
    try std.testing.expectEqual(8384, F[2].re);
    try std.testing.expectEqual(-12096, F[2].im);
    try std.testing.expectEqual(61440, F[3].re);
    try std.testing.expectEqual(-70400, F[3].im);
    try std.testing.expectEqual(245440, F[4].re);
    try std.testing.expectEqual(-369024, F[4].im);
    try std.testing.expectEqual(0, F[5].re);
    try std.testing.expectEqual(-8832, F[5].im);
    try std.testing.expectEqual(17664, F[6].re);
    try std.testing.expectEqual(-67264, F[6].im);
    try std.testing.expectEqual(423296, F[7].re);
    try std.testing.expectEqual(-556800, F[7].im);
    try std.testing.expectEqual(2350784, F[8].re);
    try std.testing.expectEqual(-2870528, F[8].im);
    try std.testing.expectEqual(12937600, F[9].re);
    try std.testing.expectEqual(-13494272, F[9].im);
    try std.testing.expectEqual(0, F[10].re);
    try std.testing.expectEqual(-96832, F[10].im);
    try std.testing.expectEqual(193664, F[11].re);
    try std.testing.expectEqual(-517120, F[11].im);
    try std.testing.expectEqual(2878144, F[12].re);
    try std.testing.expectEqual(-3966976, F[12].im);
    try std.testing.expectEqual(18542080, F[13].re);
    try std.testing.expectEqual(-20043392, F[13].im);
    try std.testing.expectEqual(50077440, F[14].re);
    try std.testing.expectEqual(-51462272, F[14].im);
    try std.testing.expectEqual(0, F[15].re);
    try std.testing.expectEqual(-478592, F[15].im);
    try std.testing.expectEqual(957184, F[16].re);
    try std.testing.expectEqual(-2589696, F[16].im);
    try std.testing.expectEqual(14495872, F[17].re);
    try std.testing.expectEqual(-17073024, F[17].im);
    try std.testing.expectEqual(52687104, F[18].re);
    try std.testing.expectEqual(-55769536, F[18].im);
    try std.testing.expectEqual(119024000, F[19].re);
    try std.testing.expectEqual(-121883392, F[19].im);

    trmm(.col_major, .right, .upper, .no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(128, F[0].re);
    try std.testing.expectEqual(-128, F[0].im);
    try std.testing.expectEqual(2560, F[1].re);
    try std.testing.expectEqual(-2048, F[1].im);
    try std.testing.expectEqual(40960, F[2].re);
    try std.testing.expectEqual(-7424, F[2].im);
    try std.testing.expectEqual(263680, F[3].re);
    try std.testing.expectEqual(-17920, F[3].im);
    try std.testing.expectEqual(1230464, F[4].re);
    try std.testing.expectEqual(-247168, F[4].im);
    try std.testing.expectEqual(45312, F[5].re);
    try std.testing.expectEqual(-14592, F[5].im);
    try std.testing.expectEqual(460160, F[6].re);
    try std.testing.expectEqual(102016, F[6].im);
    try std.testing.expectEqual(3649792, F[7].re);
    try std.testing.expectEqual(1207552, F[7].im);
    try std.testing.expectEqual(28158592, F[8].re);
    try std.testing.expectEqual(10741632, F[8].im);
    try std.testing.expectEqual(53338368, F[9].re);
    try std.testing.expectEqual(-1107712, F[9].im);
    try std.testing.expectEqual(3954560, F[10].re);
    try std.testing.expectEqual(1023104, F[10].im);
    try std.testing.expectEqual(31245568, F[11].re);
    try std.testing.expectEqual(22374656, F[11].im);
    try std.testing.expectEqual(245465984, F[12].re);
    try std.testing.expectEqual(183768704, F[12].im);
    try std.testing.expectEqual(1049432832, F[13].re);
    try std.testing.expectEqual(928512768, F[13].im);
    try std.testing.expectEqual(215399424, F[14].re);
    try std.testing.expectEqual(-1031936, F[14].im);
    try std.testing.expectEqual(80557824, F[15].re);
    try std.testing.expectEqual(45702912, F[15].im);
    try std.testing.expectEqual(684491520, F[16].re);
    try std.testing.expectEqual(510907648, F[16].im);
    try std.testing.expectEqual(3229650432, F[17].re);
    try std.testing.expectEqual(2965155328, F[17].im);
    try std.testing.expectEqual(5173135232, F[18].re);
    try std.testing.expectEqual(4803528064, F[18].im);
    try std.testing.expectEqual(630246656, F[19].re);
    try std.testing.expectEqual(54509312, F[19].im);

    trmm(.row_major, .right, .upper, .conj_no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(512, F[0].re);
    try std.testing.expectEqual(-512, F[0].im);
    try std.testing.expectEqual(72704, F[1].re);
    try std.testing.expectEqual(-58368, F[1].im);
    try std.testing.expectEqual(2213376, F[2].re);
    try std.testing.expectEqual(-453120, F[2].im);
    try std.testing.expectEqual(22427648, F[3].re);
    try std.testing.expectEqual(-1853440, F[3].im);
    try std.testing.expectEqual(146703360, F[4].re);
    try std.testing.expectEqual(-26680320, F[4].im);
    try std.testing.expectEqual(181248, F[5].re);
    try std.testing.expectEqual(-58368, F[5].im);
    try std.testing.expectEqual(13246976, F[6].re);
    try std.testing.expectEqual(2739712, F[6].im);
    try std.testing.expectEqual(205058048, F[7].re);
    try std.testing.expectEqual(65882112, F[7].im);
    try std.testing.expectEqual(2361732096, F[8].re);
    try std.testing.expectEqual(887426048, F[8].im);
    try std.testing.expectEqual(7824824320, F[9].re);
    try std.testing.expectEqual(824801280, F[9].im);
    try std.testing.expectEqual(15818240, F[10].re);
    try std.testing.expectEqual(4092416, F[10].im);
    try std.testing.expectEqual(906512384, F[11].re);
    try std.testing.expectEqual(634675200, F[11].im);
    try std.testing.expectEqual(13811544064, F[12].re);
    try std.testing.expectEqual(10284238848, F[12].im);
    try std.testing.expectEqual(94691103744, F[13].re);
    try std.testing.expectEqual(81679875072, F[13].im);
    try std.testing.expectEqual(121551441920, F[14].re);
    try std.testing.expectEqual(86119398400, F[14].im);
    try std.testing.expectEqual(322231296, F[15].re);
    try std.testing.expectEqual(182811648, F[15].im);
    try std.testing.expectEqual(19810225152, F[16].re);
    try std.testing.expectEqual(14671037440, F[16].im);
    try std.testing.expectEqual(190812244992, F[17].re);
    try std.testing.expectEqual(171085556736, F[17].im);
    try std.testing.expectEqual(599949321728, F[18].re);
    try std.testing.expectEqual(550240753152, F[18].im);
    try std.testing.expectEqual(699645327360, F[19].re);
    try std.testing.expectEqual(588992860160, F[19].im);

    trmm(.col_major, .right, .upper, .conj_no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(2048, F[0].re);
    try std.testing.expectEqual(-2048, F[0].im);
    try std.testing.expectEqual(290816, F[1].re);
    try std.testing.expectEqual(-233472, F[1].im);
    try std.testing.expectEqual(8853504, F[2].re);
    try std.testing.expectEqual(-1812480, F[2].im);
    try std.testing.expectEqual(89710592, F[3].re);
    try std.testing.expectEqual(-7413760, F[3].im);
    try std.testing.expectEqual(4107706368, F[4].re);
    try std.testing.expectEqual(-747061248, F[4].im);
    try std.testing.expectEqual(6819840, F[5].re);
    try std.testing.expectEqual(-3035136, F[5].im);
    try std.testing.expectEqual(424036352, F[6].re);
    try std.testing.expectEqual(65837056, F[6].im);
    try std.testing.expectEqual(6279888896, F[7].re);
    try std.testing.expectEqual(1800216576, F[7].im);
    try std.testing.expectEqual(129851852800, F[8].re);
    try std.testing.expectEqual(44865476608, F[8].im);
    try std.testing.expectEqual(406902763520, F[9].re);
    try std.testing.expectEqual(42884296704, F[9].im);
    try std.testing.expectEqual(1555791872, F[10].re);
    try std.testing.expectEqual(324374528, F[10].im);
    try std.testing.expectEqual(57968246784, F[11].re);
    try std.testing.expectEqual(36083900416, F[11].im);
    try std.testing.expectEqual(1229697921024, F[12].re);
    try std.testing.expectEqual(843682533376, F[12].im);
    try std.testing.expectEqual(7759928213504, F[13].re);
    try std.testing.expectEqual(6267048493056, F[13].im);
    try std.testing.expectEqual(9240090949632, F[14].re);
    try std.testing.expectEqual(6545526233088, F[14].im);
    try std.testing.expectEqual(105137786880, F[15].re);
    try std.testing.expectEqual(63951663104, F[15].im);
    try std.testing.expectEqual(3537120036864, F[16].re);
    try std.testing.expectEqual(2533685958656, F[16].im);
    try std.testing.expectEqual(28891476353024, F[17].re);
    try std.testing.expectEqual(25025695358976, F[17].im);
    try std.testing.expectEqual(71666677532672, F[18].re);
    try std.testing.expectEqual(63292117096448, F[18].im);
    try std.testing.expectEqual(70098795110400, F[19].re);
    try std.testing.expectEqual(58980867989504, F[19].im);

    trmm(.row_major, .right, .upper, .conj_no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(8192, F[0].re);
    try std.testing.expectEqual(0, F[0].im);
    try std.testing.expectEqual(1064960, F[1].re);
    try std.testing.expectEqual(98304, F[1].im);
    try std.testing.expectEqual(30662656, F[2].re);
    try std.testing.expectEqual(6586368, F[2].im);
    try std.testing.expectEqual(700547072, F[3].re);
    try std.testing.expectEqual(54657024, F[3].im);
    try std.testing.expectEqual(17429266432, F[4].re);
    try std.testing.expectEqual(6010060800, F[4].im);
    try std.testing.expectEqual(19709952, F[5].re);
    try std.testing.expectEqual(7569408, F[5].im);
    try std.testing.expectEqual(770957312, F[6].re);
    try std.testing.expectEqual(955465728, F[6].im);
    try std.testing.expectEqual(22610345984, F[7].re);
    try std.testing.expectEqual(18230575104, F[7].im);
    try std.testing.expectEqual(537020956672, F[8].re);
    try std.testing.expectEqual(452568358912, F[8].im);
    try std.testing.expectEqual(11510076342272, F[9].re);
    try std.testing.expectEqual(4599398023168, F[9].im);
    try std.testing.expectEqual(2462834688, F[10].re);
    try std.testing.expectEqual(3760332800, F[10].im);
    try std.testing.expectEqual(56215027712, F[11].re);
    try std.testing.expectEqual(190699290624, F[11].im);
    try std.testing.expectEqual(2645684174848, F[12].re);
    try std.testing.expectEqual(5305338216448, F[12].im);
    try std.testing.expectEqual(73960592572416, F[13].re);
    try std.testing.expectEqual(76604385689600, F[13].im);
    try std.testing.expectEqual(702315107483648, F[14].re);
    try std.testing.expectEqual(585005909319680, F[14].im);
    try std.testing.expectEqual(82372247552, F[15].re);
    try std.testing.expectEqual(338178899968, F[15].im);
    try std.testing.expectEqual(2847970451456, F[16].re);
    try std.testing.expectEqual(12653225295872, F[16].im);
    try std.testing.expectEqual(122181056610304, F[17].re);
    try std.testing.expectEqual(189679714058240, F[17].im);
    try std.testing.expectEqual(1763690322558976, F[18].re);
    try std.testing.expectEqual(1763592450482176, F[18].im);
    try std.testing.expectEqual(7632646195249152, F[19].re);
    try std.testing.expectEqual(6925696887062528, F[19].im);

    trmm(.col_major, .right, .upper, .conj_no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(16384, F[0].re);
    try std.testing.expectEqual(16384, F[0].im);
    try std.testing.expectEqual(1933312, F[1].re);
    try std.testing.expectEqual(2326528, F[1].im);
    try std.testing.expectEqual(48152576, F[2].re);
    try std.testing.expectEqual(74498048, F[2].im);
    try std.testing.expectEqual(1291780096, F[3].re);
    try std.testing.expectEqual(1510408192, F[3].im);
    try std.testing.expectEqual(22838607872, F[4].re);
    try std.testing.expectEqual(46878654464, F[4].im);
    try std.testing.expectEqual(49840128, F[5].re);
    try std.testing.expectEqual(56918016, F[5].im);
    try std.testing.expectEqual(366886912, F[6].re);
    try std.testing.expectEqual(3610918912, F[6].im);
    try std.testing.expectEqual(25572671488, F[7].re);
    try std.testing.expectEqual(82993610752, F[7].im);
    try std.testing.expectEqual(1005510344704, F[8].re);
    try std.testing.expectEqual(2267661549568, F[8].im);
    try std.testing.expectEqual(13822349574144, F[9].re);
    try std.testing.expectEqual(32219316387840, F[9].im);
    try std.testing.expectEqual(35760111616, F[10].re);
    try std.testing.expectEqual(58598490112, F[10].im);
    try std.testing.expectEqual(847152152576, F[11].re);
    try std.testing.expectEqual(1371301150720, F[11].im);
    try std.testing.expectEqual(34531391438848, F[12].re);
    try std.testing.expectEqual(48895650758656, F[12].im);
    try std.testing.expectEqual(823439318843392, F[13].re);
    try std.testing.expectEqual(632287135203328, F[13].im);
    try std.testing.expectEqual(234850107932672, F[14].re);
    try std.testing.expectEqual(2574978170765312, F[14].im);
    try std.testing.expectEqual(5118207229952, F[15].re);
    try std.testing.expectEqual(15814628376576, F[15].im);
    try std.testing.expectEqual(285314875244544, F[16].re);
    try std.testing.expectEqual(582480034643968, F[16].im);
    try std.testing.expectEqual(8024148419477504, F[17].re);
    try std.testing.expectEqual(8400887860035584, F[17].im);
    try std.testing.expectEqual(67422743063281660, F[18].re);
    try std.testing.expectEqual(63215563425628160, F[18].im);
    try std.testing.expectEqual(1429026691088384, F[19].re);
    try std.testing.expectEqual(29168304555556864, F[19].im);

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

    trmm(.row_major, .right, .upper, .trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-220, F[0].re);
    try std.testing.expectEqual(220, F[0].im);
    try std.testing.expectEqual(-496, F[1].re);
    try std.testing.expectEqual(496, F[1].im);
    try std.testing.expectEqual(-680, F[2].re);
    try std.testing.expectEqual(680, F[2].im);
    try std.testing.expectEqual(-704, F[3].re);
    try std.testing.expectEqual(704, F[3].im);
    try std.testing.expectEqual(-500, F[4].re);
    try std.testing.expectEqual(500, F[4].im);
    try std.testing.expectEqual(-520, F[5].re);
    try std.testing.expectEqual(520, F[5].im);
    try std.testing.expectEqual(-1176, F[6].re);
    try std.testing.expectEqual(1176, F[6].im);
    try std.testing.expectEqual(-1520, F[7].re);
    try std.testing.expectEqual(1520, F[7].im);
    try std.testing.expectEqual(-1484, F[8].re);
    try std.testing.expectEqual(1484, F[8].im);
    try std.testing.expectEqual(-1000, F[9].re);
    try std.testing.expectEqual(1000, F[9].im);
    try std.testing.expectEqual(-820, F[10].re);
    try std.testing.expectEqual(820, F[10].im);
    try std.testing.expectEqual(-1856, F[11].re);
    try std.testing.expectEqual(1856, F[11].im);
    try std.testing.expectEqual(-2360, F[12].re);
    try std.testing.expectEqual(2360, F[12].im);
    try std.testing.expectEqual(-2264, F[13].re);
    try std.testing.expectEqual(2264, F[13].im);
    try std.testing.expectEqual(-1500, F[14].re);
    try std.testing.expectEqual(1500, F[14].im);
    try std.testing.expectEqual(-1120, F[15].re);
    try std.testing.expectEqual(1120, F[15].im);
    try std.testing.expectEqual(-2536, F[16].re);
    try std.testing.expectEqual(2536, F[16].im);
    try std.testing.expectEqual(-3200, F[17].re);
    try std.testing.expectEqual(3200, F[17].im);
    try std.testing.expectEqual(-3044, F[18].re);
    try std.testing.expectEqual(3044, F[18].im);
    try std.testing.expectEqual(-2000, F[19].re);
    try std.testing.expectEqual(2000, F[19].im);

    trmm(.col_major, .right, .upper, .trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-442240, F[0].re);
    try std.testing.expectEqual(-442240, F[0].im);
    try std.testing.expectEqual(-472160, F[1].re);
    try std.testing.expectEqual(-472160, F[1].im);
    try std.testing.expectEqual(-418720, F[2].re);
    try std.testing.expectEqual(-418720, F[2].im);
    try std.testing.expectEqual(-360640, F[3].re);
    try std.testing.expectEqual(-360640, F[3].im);
    try std.testing.expectEqual(-468880, F[4].re);
    try std.testing.expectEqual(-468880, F[4].im);
    try std.testing.expectEqual(-498112, F[5].re);
    try std.testing.expectEqual(-498112, F[5].im);
    try std.testing.expectEqual(-442160, F[6].re);
    try std.testing.expectEqual(-442160, F[6].im);
    try std.testing.expectEqual(-383808, F[7].re);
    try std.testing.expectEqual(-383808, F[7].im);
    try std.testing.expectEqual(-480400, F[8].re);
    try std.testing.expectEqual(-480400, F[8].im);
    try std.testing.expectEqual(-509408, F[9].re);
    try std.testing.expectEqual(-509408, F[9].im);
    try std.testing.expectEqual(-430688, F[10].re);
    try std.testing.expectEqual(-430688, F[10].im);
    try std.testing.expectEqual(-361152, F[11].re);
    try std.testing.expectEqual(-361152, F[11].im);
    try std.testing.expectEqual(-422816, F[12].re);
    try std.testing.expectEqual(-422816, F[12].im);
    try std.testing.expectEqual(-479264, F[13].re);
    try std.testing.expectEqual(-479264, F[13].im);
    try std.testing.expectEqual(-406224, F[14].re);
    try std.testing.expectEqual(-406224, F[14].im);
    try std.testing.expectEqual(-277120, F[15].re);
    try std.testing.expectEqual(-277120, F[15].im);
    try std.testing.expectEqual(-253600, F[16].re);
    try std.testing.expectEqual(-253600, F[16].im);
    try std.testing.expectEqual(-320000, F[17].re);
    try std.testing.expectEqual(-320000, F[17].im);
    try std.testing.expectEqual(-304400, F[18].re);
    try std.testing.expectEqual(-304400, F[18].im);
    try std.testing.expectEqual(-200000, F[19].re);
    try std.testing.expectEqual(-200000, F[19].im);

    trmm(.row_major, .right, .upper, .trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(23949760, F[0].re);
    try std.testing.expectEqual(-25718720, F[0].im);
    try std.testing.expectEqual(45137280, F[1].re);
    try std.testing.expectEqual(-47025920, F[1].im);
    try std.testing.expectEqual(48328640, F[2].re);
    try std.testing.expectEqual(-50003520, F[2].im);
    try std.testing.expectEqual(37510400, F[3].re);
    try std.testing.expectEqual(-38952960, F[3].im);
    try std.testing.expectEqual(0, F[4].re);
    try std.testing.expectEqual(-1875520, F[4].im);
    try std.testing.expectEqual(26017536, F[5].re);
    try std.testing.expectEqual(-28009984, F[5].im);
    try std.testing.expectEqual(49952576, F[6].re);
    try std.testing.expectEqual(-51721216, F[6].im);
    try std.testing.expectEqual(57466880, F[7].re);
    try std.testing.expectEqual(-59002112, F[7].im);
    try std.testing.expectEqual(40752640, F[8].re);
    try std.testing.expectEqual(-42674240, F[8].im);
    try std.testing.expectEqual(0, F[9].re);
    try std.testing.expectEqual(-2037632, F[9].im);
    try std.testing.expectEqual(23755712, F[10].re);
    try std.testing.expectEqual(-25478464, F[10].im);
    try std.testing.expectEqual(47032576, F[11].re);
    try std.testing.expectEqual(-48477184, F[11].im);
    try std.testing.expectEqual(51212224, F[12].re);
    try std.testing.expectEqual(-52903488, F[12].im);
    try std.testing.expectEqual(32497920, F[13].re);
    try std.testing.expectEqual(-34414976, F[13].im);
    try std.testing.expectEqual(0, F[14].re);
    try std.testing.expectEqual(-1624896, F[14].im);
    try std.testing.expectEqual(14739200, F[15].re);
    try std.testing.expectEqual(-15847680, F[15].im);
    try std.testing.expectEqual(29198400, F[16].re);
    try std.testing.expectEqual(-30212800, F[16].im);
    try std.testing.expectEqual(29046400, F[17].re);
    try std.testing.expectEqual(-30326400, F[17].im);
    try std.testing.expectEqual(16000000, F[18].re);
    try std.testing.expectEqual(-17217600, F[18].im);
    try std.testing.expectEqual(0, F[19].re);
    try std.testing.expectEqual(-800000, F[19].im);

    trmm(.col_major, .right, .upper, .trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(7945714432, F[0].re);
    try std.testing.expectEqual(7519826176, F[0].im);
    try std.testing.expectEqual(5696197888, F[1].re);
    try std.testing.expectEqual(5140408064, F[1].im);
    try std.testing.expectEqual(4109297664, F[2].re);
    try std.testing.expectEqual(3584763392, F[2].im);
    try std.testing.expectEqual(4783425024, F[3].re);
    try std.testing.expectEqual(4389062144, F[3].im);
    try std.testing.expectEqual(8308278144, F[4].re);
    try std.testing.expectEqual(8004266112, F[4].im);
    try std.testing.expectEqual(5214802944, F[5].re);
    try std.testing.expectEqual(4761956864, F[5].im);
    try std.testing.expectEqual(3051955584, F[6].re);
    try std.testing.expectEqual(2544736896, F[6].im);
    try std.testing.expectEqual(3707885056, F[7].re);
    try std.testing.expectEqual(3256758784, F[7].im);
    try std.testing.expectEqual(6755482496, F[8].re);
    try std.testing.expectEqual(6369689728, F[8].im);
    try std.testing.expectEqual(5271982336, F[9].re);
    try std.testing.expectEqual(5008043776, F[9].im);
    try std.testing.expectEqual(1799480064, F[10].re);
    try std.testing.expectEqual(1468554496, F[10].im);
    try std.testing.expectEqual(1405652480, F[11].re);
    try std.testing.expectEqual(1058333184, F[11].im);
    try std.testing.expectEqual(3108660224, F[12].re);
    try std.testing.expectEqual(2799663872, F[12].im);
    try std.testing.expectEqual(3045160192, F[13].re);
    try std.testing.expectEqual(2784620288, F[13].im);
    try std.testing.expectEqual(1656139392, F[14].re);
    try std.testing.expectEqual(1532750208, F[14].im);
    try std.testing.expectEqual(137973760, F[15].re);
    try std.testing.expectEqual(-2216960, F[15].im);
    try std.testing.expectEqual(118822400, F[16].re);
    try std.testing.expectEqual(-2028800, F[16].im);
    try std.testing.expectEqual(118745600, F[17].re);
    try std.testing.expectEqual(-2560000, F[17].im);
    try std.testing.expectEqual(66435200, F[18].re);
    try std.testing.expectEqual(-2435200, F[18].im);
    try std.testing.expectEqual(1600000, F[19].re);
    try std.testing.expectEqual(-1600000, F[19].im);

    trmm(.row_major, .right, .upper, .conj_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(369364376064, F[0].re);
    try std.testing.expectEqual(344530046464, F[0].im);
    try std.testing.expectEqual(795525492736, F[1].re);
    try std.testing.expectEqual(736820736000, F[1].im);
    try std.testing.expectEqual(980051968512, F[2].re);
    try std.testing.expectEqual(912451143168, F[2].im);
    try std.testing.expectEqual(1028202553344, F[3].re);
    try std.testing.expectEqual(973910011904, F[3].im);
    try std.testing.expectEqual(830827814400, F[4].re);
    try std.testing.expectEqual(800426611200, F[4].im);
    try std.testing.expectEqual(303296843776, F[5].re);
    try std.testing.expectEqual(280562739200, F[5].im);
    try std.testing.expectEqual(658183741440, F[6].re);
    try std.testing.expectEqual(605099495424, F[6].im);
    try std.testing.expectEqual(887435982848, F[7].re);
    try std.testing.expectEqual(826536708096, F[7].im);
    try std.testing.expectEqual(935175256576, F[8].re);
    try std.testing.expectEqual(884739921408, F[8].im);
    try std.testing.expectEqual(527198233600, F[9].re);
    try std.testing.expectEqual(500804377600, F[9].im);
    try std.testing.expectEqual(137592413696, F[10].re);
    try std.testing.expectEqual(123145778688, F[10].im);
    try std.testing.expectEqual(314706739200, F[11].re);
    try std.testing.expectEqual(280778911744, F[11].im);
    try std.testing.expectEqual(431547665920, F[12].re);
    try std.testing.expectEqual(393486269952, F[12].im);
    try std.testing.expectEqual(363923325952, F[13].re);
    try std.testing.expectEqual(334251158528, F[13].im);
    try std.testing.expectEqual(165613939200, F[14].re);
    try std.testing.expectEqual(153275020800, F[14].im);
    try std.testing.expectEqual(4022384640, F[15].re);
    try std.testing.expectEqual(-126781440, F[15].im);
    try std.testing.expectEqual(9582553600, F[16].re);
    try std.testing.expectEqual(-290393600, F[16].im);
    try std.testing.expectEqual(9991142400, F[17].re);
    try std.testing.expectEqual(-365491200, F[17].im);
    try std.testing.expectEqual(5177075200, F[18].re);
    try std.testing.expectEqual(-313075200, F[18].im);
    try std.testing.expectEqual(160000000, F[19].re);
    try std.testing.expectEqual(-160000000, F[19].im);

    trmm(.col_major, .right, .upper, .conj_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(90989021460480, F[0].re);
    try std.testing.expectEqual(84675643611136, F[0].im);
    try std.testing.expectEqual(57788297322496, F[1].re);
    try std.testing.expectEqual(53077554184192, F[1].im);
    try std.testing.expectEqual(36804850296832, F[2].re);
    try std.testing.expectEqual(33373909739520, F[2].im);
    try std.testing.expectEqual(39529242943488, F[3].re);
    try std.testing.expectEqual(36065239146496, F[3].im);
    try std.testing.expectEqual(98340097118208, F[4].re);
    try std.testing.expectEqual(91610973061120, F[4].im);
    try std.testing.expectEqual(59423833534464, F[5].re);
    try std.testing.expectEqual(54591282376704, F[5].im);
    try std.testing.expectEqual(36750911100928, F[6].re);
    try std.testing.expectEqual(33248934045696, F[6].im);
    try std.testing.expectEqual(40241733156864, F[7].re);
    try std.testing.expectEqual(36597714452480, F[7].im);
    try std.testing.expectEqual(80582140219392, F[8].re);
    try std.testing.expectEqual(74310771138560, F[8].im);
    try std.testing.expectEqual(54535972716544, F[9].re);
    try std.testing.expectEqual(50074285858816, F[9].im);
    try std.testing.expectEqual(19555300052992, F[10].re);
    try std.testing.expectEqual(17410579070976, F[10].im);
    try std.testing.expectEqual(16669082132480, F[11].re);
    try std.testing.expectEqual(14576655147008, F[11].im);
    try std.testing.expectEqual(33717547755520, F[12].re);
    try std.testing.expectEqual(29877078730752, F[12].im);
    try std.testing.expectEqual(28617322442752, F[13].re);
    try std.testing.expectEqual(25368000892928, F[13].im);
    try std.testing.expectEqual(13083658598400, F[14].re);
    try std.testing.expectEqual(11618846361600, F[14].im);
    try std.testing.expectEqual(321061232640, F[15].re);
    try std.testing.expectEqual(-24995389440, F[15].im);
    try std.testing.expectEqual(958255360000, F[16].re);
    try std.testing.expectEqual(-29039360000, F[16].im);
    try std.testing.expectEqual(999114240000, F[17].re);
    try std.testing.expectEqual(-36549120000, F[17].im);
    try std.testing.expectEqual(517707520000, F[18].re);
    try std.testing.expectEqual(-31307520000, F[18].im);
    try std.testing.expectEqual(16000000000, F[19].re);
    try std.testing.expectEqual(-16000000000, F[19].im);

    trmm(.row_major, .right, .upper, .conj_trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(3515861167300608, F[0].re);
    try std.testing.expectEqual(3585699968057344, F[0].im);
    try std.testing.expectEqual(6543833326469120, F[1].re);
    try std.testing.expectEqual(6252484346396672, F[1].im);
    try std.testing.expectEqual(8120905313042432, F[2].re);
    try std.testing.expectEqual(7656669295943680, F[2].im);
    try std.testing.expectEqual(7874135777050624, F[3].re);
    try std.testing.expectEqual(7480066809069568, F[3].im);
    try std.testing.expectEqual(13458248114176, F[4].re);
    try std.testing.expectEqual(379902140358656, F[4].im);
    try std.testing.expectEqual(3166606886846464, F[5].re);
    try std.testing.expectEqual(3123652333010944, F[5].im);
    try std.testing.expectEqual(6377135371689984, F[6].re);
    try std.testing.expectEqual(5989285748113408, F[6].im);
    try std.testing.expectEqual(7792046252687360, F[7].re);
    try std.testing.expectEqual(7319539230507008, F[7].im);
    try std.testing.expectEqual(4375420555485184, F[8].re);
    try std.testing.expectEqual(4315728691421184, F[8].im);
    try std.testing.expectEqual(8923373715456, F[9].re);
    try std.testing.expectEqual(209220517150720, F[9].im);
    try std.testing.expectEqual(1261803003142144, F[10].re);
    try std.testing.expectEqual(1187334885711872, F[10].im);
    try std.testing.expectEqual(2636716334022656, F[11].re);
    try std.testing.expectEqual(2396559880552448, F[11].im);
    try std.testing.expectEqual(2395270510747648, F[12].re);
    try std.testing.expectEqual(2244928084672512, F[12].im);
    try std.testing.expectEqual(1053191330971648, F[13].re);
    try std.testing.expectEqual(1037478355599360, F[13].im);
    try std.testing.expectEqual(2929624473600, F[14].re);
    try std.testing.expectEqual(49405009920000, F[14].im);
    try std.testing.expectEqual(28950847324160, F[15].re);
    try std.testing.expectEqual(-899692953600, F[15].im);
    try std.testing.expectEqual(53223715840000, F[16].re);
    try std.testing.expectEqual(-1078210560000, F[16].im);
    try std.testing.expectEqual(32022947840000, F[17].re);
    try std.testing.expectEqual(-788090880000, F[17].im);
    try std.testing.expectEqual(2378030080000, F[18].re);
    try std.testing.expectEqual(-307200000000, F[18].im);
    try std.testing.expectEqual(64000000000, F[19].re);
    try std.testing.expectEqual(0, F[19].im);

    trmm(.col_major, .right, .upper, .conj_trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(350469929612984300, F[0].re);
    try std.testing.expectEqual(356797663793856500, F[0].im);
    try std.testing.expectEqual(147068064488685570, F[1].re);
    try std.testing.expectEqual(176098409217064960, F[1].im);
    try std.testing.expectEqual(209886303586041860, F[2].re);
    try std.testing.expectEqual(230676857978896400, F[2].im);
    try std.testing.expectEqual(305670996926201860, F[3].re);
    try std.testing.expectEqual(311768401099685900, F[3].im);
    try std.testing.expectEqual(376849380603559940, F[4].re);
    try std.testing.expectEqual(360501925193613300, F[4].im);
    try std.testing.expectEqual(74949260962005000, F[5].re);
    try std.testing.expectEqual(93102279446265860, F[5].im);
    try std.testing.expectEqual(61750724509220860, F[6].re);
    try std.testing.expectEqual(85057423828336640, F[6].im);
    try std.testing.expectEqual(129481687695491070, F[7].re);
    try std.testing.expectEqual(145196866112061440, F[7].im);
    try std.testing.expectEqual(177475442359238660, F[8].re);
    try std.testing.expectEqual(178917925218713600, F[8].im);
    try std.testing.expectEqual(78375292744368130, F[9].re);
    try std.testing.expectEqual(75062225023926270, F[9].im);
    try std.testing.expectEqual(578647964319744, F[10].re);
    try std.testing.expectEqual(8427174091948032, F[10].im);
    try std.testing.expectEqual(2570661914279936, F[11].re);
    try std.testing.expectEqual(10001774536491008, F[11].im);
    try std.testing.expectEqual(5410161572790272, F[12].re);
    try std.testing.expectEqual(9176888977080320, F[12].im);
    try std.testing.expectEqual(3105628943384576, F[13].re);
    try std.testing.expectEqual(4105682648662016, F[13].im);
    try std.testing.expectEqual(135340116787200, F[14].re);
    try std.testing.expectEqual(75178068787200, F[14].im);
    try std.testing.expectEqual(65845080555520, F[15].re);
    try std.testing.expectEqual(56102308741120, F[15].im);
    try std.testing.expectEqual(108603852800000, F[16].re);
    try std.testing.expectEqual(104291010560000, F[16].im);
    try std.testing.expectEqual(65622077440000, F[17].re);
    try std.testing.expectEqual(62469713920000, F[17].im);
    try std.testing.expectEqual(5370460160000, F[18].re);
    try std.testing.expectEqual(4141660160000, F[18].im);
    try std.testing.expectEqual(128000000000, F[19].re);
    try std.testing.expectEqual(128000000000, F[19].im);

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

    trmm(.row_major, .right, .lower, .no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-860, F[0].re);
    try std.testing.expectEqual(860, F[0].im);
    try std.testing.expectEqual(-912, F[1].re);
    try std.testing.expectEqual(912, F[1].im);
    try std.testing.expectEqual(-904, F[2].re);
    try std.testing.expectEqual(904, F[2].im);
    try std.testing.expectEqual(-784, F[3].re);
    try std.testing.expectEqual(784, F[3].im);
    try std.testing.expectEqual(-500, F[4].re);
    try std.testing.expectEqual(500, F[4].im);
    try std.testing.expectEqual(-1960, F[5].re);
    try std.testing.expectEqual(1960, F[5].im);
    try std.testing.expectEqual(-2072, F[6].re);
    try std.testing.expectEqual(2072, F[6].im);
    try std.testing.expectEqual(-1984, F[7].re);
    try std.testing.expectEqual(1984, F[7].im);
    try std.testing.expectEqual(-1644, F[8].re);
    try std.testing.expectEqual(1644, F[8].im);
    try std.testing.expectEqual(-1000, F[9].re);
    try std.testing.expectEqual(1000, F[9].im);
    try std.testing.expectEqual(-3060, F[10].re);
    try std.testing.expectEqual(3060, F[10].im);
    try std.testing.expectEqual(-3232, F[11].re);
    try std.testing.expectEqual(3232, F[11].im);
    try std.testing.expectEqual(-3064, F[12].re);
    try std.testing.expectEqual(3064, F[12].im);
    try std.testing.expectEqual(-2504, F[13].re);
    try std.testing.expectEqual(2504, F[13].im);
    try std.testing.expectEqual(-1500, F[14].re);
    try std.testing.expectEqual(1500, F[14].im);
    try std.testing.expectEqual(-4160, F[15].re);
    try std.testing.expectEqual(4160, F[15].im);
    try std.testing.expectEqual(-4392, F[16].re);
    try std.testing.expectEqual(4392, F[16].im);
    try std.testing.expectEqual(-4144, F[17].re);
    try std.testing.expectEqual(4144, F[17].im);
    try std.testing.expectEqual(-3364, F[18].re);
    try std.testing.expectEqual(3364, F[18].im);
    try std.testing.expectEqual(-2000, F[19].re);
    try std.testing.expectEqual(2000, F[19].im);

    trmm(.col_major, .right, .lower, .no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-164032, F[0].re);
    try std.testing.expectEqual(-164032, F[0].im);
    try std.testing.expectEqual(-154272, F[1].re);
    try std.testing.expectEqual(-154272, F[1].im);
    try std.testing.expectEqual(-148192, F[2].re);
    try std.testing.expectEqual(-148192, F[2].im);
    try std.testing.expectEqual(-164352, F[3].re);
    try std.testing.expectEqual(-164352, F[3].im);
    try std.testing.expectEqual(-352592, F[4].re);
    try std.testing.expectEqual(-352592, F[4].im);
    try std.testing.expectEqual(-342784, F[5].re);
    try std.testing.expectEqual(-342784, F[5].im);
    try std.testing.expectEqual(-344496, F[6].re);
    try std.testing.expectEqual(-344496, F[6].im);
    try std.testing.expectEqual(-388736, F[7].re);
    try std.testing.expectEqual(-388736, F[7].im);
    try std.testing.expectEqual(-520592, F[8].re);
    try std.testing.expectEqual(-520592, F[8].im);
    try std.testing.expectEqual(-440864, F[9].re);
    try std.testing.expectEqual(-440864, F[9].im);
    try std.testing.expectEqual(-444960, F[10].re);
    try std.testing.expectEqual(-444960, F[10].im);
    try std.testing.expectEqual(-521024, F[11].re);
    try std.testing.expectEqual(-521024, F[11].im);
    try std.testing.expectEqual(-584224, F[12].re);
    try std.testing.expectEqual(-584224, F[12].im);
    try std.testing.expectEqual(-521824, F[13].re);
    try std.testing.expectEqual(-521824, F[13].im);
    try std.testing.expectEqual(-383120, F[14].re);
    try std.testing.expectEqual(-383120, F[14].im);
    try std.testing.expectEqual(-476160, F[15].re);
    try std.testing.expectEqual(-476160, F[15].im);
    try std.testing.expectEqual(-439200, F[16].re);
    try std.testing.expectEqual(-439200, F[16].im);
    try std.testing.expectEqual(-414400, F[17].re);
    try std.testing.expectEqual(-414400, F[17].im);
    try std.testing.expectEqual(-336400, F[18].re);
    try std.testing.expectEqual(-336400, F[18].im);
    try std.testing.expectEqual(-200000, F[19].re);
    try std.testing.expectEqual(-200000, F[19].im);

    trmm(.row_major, .right, .lower, .no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(50359232, F[0].re);
    try std.testing.expectEqual(-51015360, F[0].im);
    try std.testing.expectEqual(49317248, F[1].re);
    try std.testing.expectEqual(-49934336, F[1].im);
    try std.testing.expectEqual(44271808, F[2].re);
    try std.testing.expectEqual(-44864576, F[2].im);
    try std.testing.expectEqual(33848832, F[3].re);
    try std.testing.expectEqual(-34506240, F[3].im);
    try std.testing.expectEqual(0, F[4].re);
    try std.testing.expectEqual(-1410368, F[4].im);
    try std.testing.expectEqual(95722752, F[5].re);
    try std.testing.expectEqual(-97093888, F[5].im);
    try std.testing.expectEqual(92855616, F[6].re);
    try std.testing.expectEqual(-94233600, F[6].im);
    try std.testing.expectEqual(78042112, F[7].re);
    try std.testing.expectEqual(-79597056, F[7].im);
    try std.testing.expectEqual(42322944, F[8].re);
    try std.testing.expectEqual(-44405312, F[8].im);
    try std.testing.expectEqual(0, F[9].re);
    try std.testing.expectEqual(-1763456, F[9].im);
    try std.testing.expectEqual(103789248, F[10].re);
    try std.testing.expectEqual(-105569088, F[10].im);
    try std.testing.expectEqual(97241344, F[11].re);
    try std.testing.expectEqual(-99325440, F[11].im);
    try std.testing.expectEqual(72818368, F[12].re);
    try std.testing.expectEqual(-75155264, F[12].im);
    try std.testing.expectEqual(36779520, F[13].re);
    try std.testing.expectEqual(-38866816, F[13].im);
    try std.testing.expectEqual(0, F[14].re);
    try std.testing.expectEqual(-1532480, F[14].im);
    try std.testing.expectEqual(67104000, F[15].re);
    try std.testing.expectEqual(-69008640, F[15].im);
    try std.testing.expectEqual(60366400, F[16].re);
    try std.testing.expectEqual(-62123200, F[16].im);
    try std.testing.expectEqual(42620800, F[17].re);
    try std.testing.expectEqual(-44278400, F[17].im);
    try std.testing.expectEqual(19200000, F[18].re);
    try std.testing.expectEqual(-20545600, F[18].im);
    try std.testing.expectEqual(0, F[19].re);
    try std.testing.expectEqual(-800000, F[19].im);

    trmm(.col_major, .right, .lower, .no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(3191844096, F[0].re);
    try std.testing.expectEqual(2878984960, F[0].im);
    try std.testing.expectEqual(2503852800, F[1].re);
    try std.testing.expectEqual(2205436160, F[1].im);
    try std.testing.expectEqual(2634402304, F[2].re);
    try std.testing.expectEqual(2371130368, F[2].im);
    try std.testing.expectEqual(3085530112, F[3].re);
    try std.testing.expectEqual(2863582208, F[3].im);
    try std.testing.expectEqual(6614308224, F[4].re);
    try std.testing.expectEqual(6387630720, F[4].im);
    try std.testing.expectEqual(3612405248, F[5].re);
    try std.testing.expectEqual(3026152448, F[5].im);
    try std.testing.expectEqual(4629382528, F[6].re);
    try std.testing.expectEqual(4086499968, F[6].im);
    try std.testing.expectEqual(6010003456, F[7].re);
    try std.testing.expectEqual(5524357120, F[7].im);
    try std.testing.expectEqual(8109543296, F[8].re);
    try std.testing.expectEqual(7695647872, F[8].im);
    try std.testing.expectEqual(4836772608, F[9].re);
    try std.testing.expectEqual(4613374208, F[9].im);
    try std.testing.expectEqual(1737271552, F[10].re);
    try std.testing.expectEqual(1148440320, F[10].im);
    try std.testing.expectEqual(4305617408, F[11].re);
    try std.testing.expectEqual(3753655808, F[11].im);
    try std.testing.expectEqual(5265803264, F[12].re);
    try std.testing.expectEqual(4824638208, F[12].im);
    try std.testing.expectEqual(3693564672, F[13].re);
    try std.testing.expectEqual(3405489408, F[13].im);
    try std.testing.expectEqual(1646712960, F[14].re);
    try std.testing.expectEqual(1532935040, F[14].im);
    try std.testing.expectEqual(336225280, F[15].re);
    try std.testing.expectEqual(-3809280, F[15].im);
    try std.testing.expectEqual(244979200, F[16].re);
    try std.testing.expectEqual(-3513600, F[16].im);
    try std.testing.expectEqual(173798400, F[17].re);
    try std.testing.expectEqual(-3315200, F[17].im);
    try std.testing.expectEqual(79491200, F[18].re);
    try std.testing.expectEqual(-2691200, F[18].im);
    try std.testing.expectEqual(1600000, F[19].re);
    try std.testing.expectEqual(-1600000, F[19].im);

    trmm(.row_major, .right, .lower, .conj_no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(941849362944, F[0].re);
    try std.testing.expectEqual(888606385664, F[0].im);
    try std.testing.expectEqual(988434360320, F[1].re);
    try std.testing.expectEqual(932401563648, F[1].im);
    try std.testing.expectEqual(967663444480, F[2].re);
    try std.testing.expectEqual(917138724352, F[2].im);
    try std.testing.expectEqual(869473878016, F[3].re);
    try std.testing.expectEqual(830844796928, F[3].im);
    try std.testing.expectEqual(661430822400, F[4].re);
    try std.testing.expectEqual(638763072000, F[4].im);
    try std.testing.expectEqual(1315294623744, F[5].re);
    try std.testing.expectEqual(1233297219584, F[5].im);
    try std.testing.expectEqual(1395187810304, F[6].re);
    try std.testing.expectEqual(1308872126464, F[6].im);
    try std.testing.expectEqual(1341390376960, F[7].re);
    try std.testing.expectEqual(1265783644160, F[7].im);
    try std.testing.expectEqual(1080655460864, F[8].re);
    try std.testing.expectEqual(1027753162240, F[8].im);
    try std.testing.expectEqual(483677260800, F[9].re);
    try std.testing.expectEqual(461337420800, F[9].im);
    try std.testing.expectEqual(716691275264, F[10].re);
    try std.testing.expectEqual(653683447296, F[10].im);
    try std.testing.expectEqual(769388982272, F[11].re);
    try std.testing.expectEqual(703156559872, F[11].im);
    try std.testing.expectEqual(691256018432, F[12].re);
    try std.testing.expectEqual(637106447872, F[12].im);
    try std.testing.expectEqual(438795359232, F[13].re);
    try std.testing.expectEqual(405978958848, F[13].im);
    try std.testing.expectEqual(164671296000, F[14].re);
    try std.testing.expectEqual(153293504000, F[14].im);
    try std.testing.expectEqual(20093368320, F[15].re);
    try std.testing.expectEqual(-552069120, F[15].im);
    try std.testing.expectEqual(20747942400, F[16].re);
    try std.testing.expectEqual(-581312000, F[16].im);
    try std.testing.expectEqual(14908083200, F[17].re);
    try std.testing.expectEqual(-513356800, F[17].im);
    try std.testing.expectEqual(6194931200, F[18].re);
    try std.testing.expectEqual(-358131200, F[18].im);
    try std.testing.expectEqual(160000000, F[19].re);
    try std.testing.expectEqual(-160000000, F[19].im);

    trmm(.col_major, .right, .lower, .conj_no_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(33501764704256, F[0].re);
    try std.testing.expectEqual(31179644991488, F[0].im);
    try std.testing.expectEqual(27599108972544, F[1].re);
    try std.testing.expectEqual(25617429266432, F[1].im);
    try std.testing.expectEqual(26391090923520, F[2].re);
    try std.testing.expectEqual(24429266716672, F[2].im);
    try std.testing.expectEqual(23766380208128, F[3].re);
    try std.testing.expectEqual(21875493953536, F[3].im);
    try std.testing.expectEqual(78816172134400, F[4].re);
    try std.testing.expectEqual(73686046851072, F[4].im);
    try std.testing.expectEqual(68698878070784, F[5].re);
    try std.testing.expectEqual(63889827860480, F[5].im);
    try std.testing.expectEqual(68175343400960, F[6].re);
    try std.testing.expectEqual(63070530750464, F[6].im);
    try std.testing.expectEqual(62909139247104, F[7].re);
    try std.testing.expectEqual(57916677464064, F[7].im);
    try std.testing.expectEqual(96149297541120, F[8].re);
    try std.testing.expectEqual(89086246797312, F[8].im);
    try std.testing.expectEqual(50618242670592, F[9].re);
    try std.testing.expectEqual(46693566169088, F[9].im);
    try std.testing.expectEqual(46861234761728, F[10].re);
    try std.testing.expectEqual(42554487611392, F[10].im);
    try std.testing.expectEqual(41143055704064, F[11].re);
    try std.testing.expectEqual(36523625242624, F[11].im);
    try std.testing.expectEqual(54195292792832, F[12].re);
    try std.testing.expectEqual(48373585078272, F[12].im);
    try std.testing.expectEqual(34541093957632, F[13].re);
    try std.testing.expectEqual(30813332328448, F[13].im);
    try std.testing.expectEqual(13010612992000, F[14].re);
    try std.testing.expectEqual(11621655808000, F[14].im);
    try std.testing.expectEqual(1539895992320, F[15].re);
    try std.testing.expectEqual(-54757253120, F[15].im);
    try std.testing.expectEqual(2074794240000, F[16].re);
    try std.testing.expectEqual(-58131200000, F[16].im);
    try std.testing.expectEqual(1490808320000, F[17].re);
    try std.testing.expectEqual(-51335680000, F[17].im);
    try std.testing.expectEqual(619493120000, F[18].re);
    try std.testing.expectEqual(-35813120000, F[18].im);
    try std.testing.expectEqual(16000000000, F[19].re);
    try std.testing.expectEqual(-16000000000, F[19].im);

    trmm(.row_major, .right, .lower, .conj_no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(9969837648011264, F[0].re);
    try std.testing.expectEqual(9408728405835776, F[0].im);
    try std.testing.expectEqual(9822672725721088, F[1].re);
    try std.testing.expectEqual(9250943590612992, F[1].im);
    try std.testing.expectEqual(8966190859763712, F[2].re);
    try std.testing.expectEqual(8455792590233600, F[2].im);
    try std.testing.expectEqual(7570134297411584, F[3].re);
    try std.testing.expectEqual(7165144246026240, F[3].im);
    try std.testing.expectEqual(10260250566656, F[4].re);
    try std.testing.expectEqual(305004437970944, F[4].im);
    try std.testing.expectEqual(14819315895877632, F[5].re);
    try std.testing.expectEqual(13950983311523840, F[5].im);
    try std.testing.expectEqual(14022405896970240, F[6].re);
    try std.testing.expectEqual(13209390871674880, F[6].im);
    try std.testing.expectEqual(11589612672221184, F[7].re);
    try std.testing.expectEqual(10951669490384896, F[7].im);
    try std.testing.expectEqual(4873477397864448, F[8].re);
    try std.testing.expectEqual(4853053440909312, F[8].im);
    try std.testing.expectEqual(7849353003008, F[9].re);
    try std.testing.expectEqual(194623617679360, F[9].im);
    try std.testing.expectEqual(6684161218699264, F[10].re);
    try std.testing.expectEqual(6132108550905856, F[10].im);
    try std.testing.expectEqual(6104341247393792, F[11].re);
    try std.testing.expectEqual(5595277755088896, F[11].im);
    try std.testing.expectEqual(3695578575642624, F[12].re);
    try std.testing.expectEqual(3492890017726464, F[12].im);
    try std.testing.expectEqual(1256474370490368, F[13].re);
    try std.testing.expectEqual(1246387810140160, F[13].im);
    try std.testing.expectEqual(2777914368000, F[14].re);
    try std.testing.expectEqual(49264537600000, F[14].im);
    try std.testing.expectEqual(159571494010880, F[15].re);
    try std.testing.expectEqual(-4319680921600, F[15].im);
    try std.testing.expectEqual(119358182400000, F[16].re);
    try std.testing.expectEqual(-2274078720000, F[16].im);
    try std.testing.expectEqual(49159792640000, F[17].re);
    try std.testing.expectEqual(-1171599360000, F[17].im);
    try std.testing.expectEqual(2846612480000, F[18].re);
    try std.testing.expectEqual(-368640000000, F[18].im);
    try std.testing.expectEqual(64000000000, F[19].re);
    try std.testing.expectEqual(0, F[19].im);

    trmm(.col_major, .right, .lower, .conj_no_trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(121202450121539580, F[0].re);
    try std.testing.expectEqual(155274567611596800, F[0].im);
    try std.testing.expectEqual(140878963453919230, F[1].re);
    try std.testing.expectEqual(172009355512053760, F[1].im);
    try std.testing.expectEqual(193511357218701300, F[2].re);
    try std.testing.expectEqual(214885256285863940, F[2].im);
    try std.testing.expectEqual(169333400353439740, F[3].re);
    try std.testing.expectEqual(184158131176275970, F[3].im);
    try std.testing.expectEqual(293176944375988200, F[4].re);
    try std.testing.expectEqual(281581316975525900, F[4].im);
    try std.testing.expectEqual(49187313508057090, F[5].re);
    try std.testing.expectEqual(108591651371188220, F[5].im);
    try std.testing.expectEqual(215733058465415170, F[6].re);
    try std.testing.expectEqual(252449844919877630, F[6].im);
    try std.testing.expectEqual(202361940064665600, F[7].re);
    try std.testing.expectEqual(223975943974879230, F[7].im);
    try std.testing.expectEqual(214154739093897200, F[8].re);
    try std.testing.expectEqual(214918457947029500, F[8].im);
    try std.testing.expectEqual(72938603776507900, F[9].re);
    try std.testing.expectEqual(70132367347613700, F[9].im);
    try std.testing.expectEqual(1430465288994816, F[10].re);
    try std.testing.expectEqual(28369235244810240, F[10].im);
    try std.testing.expectEqual(9957970649219072, F[11].re);
    try std.testing.expectEqual(23157335873355776, F[11].im);
    try std.testing.expectEqual(9954031707832320, F[12].re);
    try std.testing.expectEqual(14195010889138176, F[12].im);
    try std.testing.expectEqual(3952956531900416, F[13].re);
    try std.testing.expectEqual(4911996412461056, F[13].im);
    try std.testing.expectEqual(134755751936000, F[14].re);
    try std.testing.expectEqual(74593703936000, F[14].im);
    try std.testing.expectEqual(332902349864960, F[15].re);
    try std.testing.expectEqual(310503626178560, F[15].im);
    try std.testing.expectEqual(243264522240000, F[16].re);
    try std.testing.expectEqual(234168207360000, F[16].im);
    try std.testing.expectEqual(100662784000000, F[17].re);
    try std.testing.expectEqual(95976386560000, F[17].im);
    try std.testing.expectEqual(6430504960000, F[18].re);
    try std.testing.expectEqual(4955944960000, F[18].im);
    try std.testing.expectEqual(128000000000, F[19].re);
    try std.testing.expectEqual(128000000000, F[19].im);

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

    trmm(.row_major, .right, .lower, .trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-4, F[0].re);
    try std.testing.expectEqual(4, F[0].im);
    try std.testing.expectEqual(-80, F[1].re);
    try std.testing.expectEqual(80, F[1].im);
    try std.testing.expectEqual(-296, F[2].re);
    try std.testing.expectEqual(296, F[2].im);
    try std.testing.expectEqual(-720, F[3].re);
    try std.testing.expectEqual(720, F[3].im);
    try std.testing.expectEqual(-1420, F[4].re);
    try std.testing.expectEqual(1420, F[4].im);
    try std.testing.expectEqual(-24, F[5].re);
    try std.testing.expectEqual(24, F[5].im);
    try std.testing.expectEqual(-340, F[6].re);
    try std.testing.expectEqual(340, F[6].im);
    try std.testing.expectEqual(-1016, F[7].re);
    try std.testing.expectEqual(1016, F[7].im);
    try std.testing.expectEqual(-2120, F[8].re);
    try std.testing.expectEqual(2120, F[8].im);
    try std.testing.expectEqual(-3720, F[9].re);
    try std.testing.expectEqual(3720, F[9].im);
    try std.testing.expectEqual(-44, F[10].re);
    try std.testing.expectEqual(44, F[10].im);
    try std.testing.expectEqual(-600, F[11].re);
    try std.testing.expectEqual(600, F[11].im);
    try std.testing.expectEqual(-1736, F[12].re);
    try std.testing.expectEqual(1736, F[12].im);
    try std.testing.expectEqual(-3520, F[13].re);
    try std.testing.expectEqual(3520, F[13].im);
    try std.testing.expectEqual(-6020, F[14].re);
    try std.testing.expectEqual(6020, F[14].im);
    try std.testing.expectEqual(-64, F[15].re);
    try std.testing.expectEqual(64, F[15].im);
    try std.testing.expectEqual(-860, F[16].re);
    try std.testing.expectEqual(860, F[16].im);
    try std.testing.expectEqual(-2456, F[17].re);
    try std.testing.expectEqual(2456, F[17].im);
    try std.testing.expectEqual(-4920, F[18].re);
    try std.testing.expectEqual(4920, F[18].im);
    try std.testing.expectEqual(-8320, F[19].re);
    try std.testing.expectEqual(8320, F[19].im);

    trmm(.col_major, .right, .lower, .trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-16, F[0].re);
    try std.testing.expectEqual(-16, F[0].im);
    try std.testing.expectEqual(-320, F[1].re);
    try std.testing.expectEqual(-320, F[1].im);
    try std.testing.expectEqual(-1184, F[2].re);
    try std.testing.expectEqual(-1184, F[2].im);
    try std.testing.expectEqual(-2880, F[3].re);
    try std.testing.expectEqual(-2880, F[3].im);
    try std.testing.expectEqual(-39792, F[4].re);
    try std.testing.expectEqual(-39792, F[4].im);
    try std.testing.expectEqual(-1312, F[5].re);
    try std.testing.expectEqual(-1312, F[5].im);
    try std.testing.expectEqual(-11888, F[6].re);
    try std.testing.expectEqual(-11888, F[6].im);
    try std.testing.expectEqual(-34208, F[7].re);
    try std.testing.expectEqual(-34208, F[7].im);
    try std.testing.expectEqual(-155728, F[8].re);
    try std.testing.expectEqual(-155728, F[8].im);
    try std.testing.expectEqual(-195168, F[9].re);
    try std.testing.expectEqual(-195168, F[9].im);
    try std.testing.expectEqual(-16720, F[10].re);
    try std.testing.expectEqual(-16720, F[10].im);
    try std.testing.expectEqual(-72352, F[11].re);
    try std.testing.expectEqual(-72352, F[11].im);
    try std.testing.expectEqual(-301840, F[12].re);
    try std.testing.expectEqual(-301840, F[12].im);
    try std.testing.expectEqual(-477984, F[13].re);
    try std.testing.expectEqual(-477984, F[13].im);
    try std.testing.expectEqual(-476960, F[14].re);
    try std.testing.expectEqual(-476960, F[14].im);
    try std.testing.expectEqual(-86560, F[15].re);
    try std.testing.expectEqual(-86560, F[15].im);
    try std.testing.expectEqual(-408960, F[16].re);
    try std.testing.expectEqual(-408960, F[16].im);
    try std.testing.expectEqual(-752960, F[17].re);
    try std.testing.expectEqual(-752960, F[17].im);
    try std.testing.expectEqual(-995760, F[18].re);
    try std.testing.expectEqual(-995760, F[18].im);
    try std.testing.expectEqual(-928160, F[19].re);
    try std.testing.expectEqual(-928160, F[19].im);

    trmm(.row_major, .right, .lower, .trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(0, F[0].re);
    try std.testing.expectEqual(-64, F[0].im);
    try std.testing.expectEqual(384, F[1].re);
    try std.testing.expectEqual(-1664, F[1].im);
    try std.testing.expectEqual(16064, F[2].re);
    try std.testing.expectEqual(-20800, F[2].im);
    try std.testing.expectEqual(108032, F[3].re);
    try std.testing.expectEqual(-119552, F[3].im);
    try std.testing.expectEqual(414912, F[4].re);
    try std.testing.expectEqual(-574080, F[4].im);
    try std.testing.expectEqual(0, F[5].re);
    try std.testing.expectEqual(-5248, F[5].im);
    try std.testing.expectEqual(31488, F[6].re);
    try std.testing.expectEqual(-79040, F[6].im);
    try std.testing.expectEqual(628352, F[7].re);
    try std.testing.expectEqual(-765184, F[7].im);
    try std.testing.expectEqual(3355328, F[8].re);
    try std.testing.expectEqual(-3978240, F[8].im);
    try std.testing.expectEqual(19253376, F[9].re);
    try std.testing.expectEqual(-20034048, F[9].im);
    try std.testing.expectEqual(0, F[10].re);
    try std.testing.expectEqual(-66880, F[10].im);
    try std.testing.expectEqual(401280, F[11].re);
    try std.testing.expectEqual(-690688, F[11].im);
    try std.testing.expectEqual(4208576, F[12].re);
    try std.testing.expectEqual(-5415936, F[12].im);
    try std.testing.expectEqual(27722496, F[13].re);
    try std.testing.expectEqual(-29634432, F[13].im);
    try std.testing.expectEqual(81427200, F[14].re);
    try std.testing.expectEqual(-83335040, F[14].im);
    try std.testing.expectEqual(0, F[15].re);
    try std.testing.expectEqual(-346240, F[15].im);
    try std.testing.expectEqual(2077440, F[16].re);
    try std.testing.expectEqual(-3713280, F[16].im);
    try std.testing.expectEqual(23438720, F[17].re);
    try std.testing.expectEqual(-26450560, F[17].im);
    try std.testing.expectEqual(87562240, F[18].re);
    try std.testing.expectEqual(-91545280, F[18].im);
    try std.testing.expectEqual(208124800, F[19].re);
    try std.testing.expectEqual(-211837440, F[19].im);

    trmm(.col_major, .right, .lower, .trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(128, F[0].re);
    try std.testing.expectEqual(-128, F[0].im);
    try std.testing.expectEqual(4096, F[1].re);
    try std.testing.expectEqual(-2560, F[1].im);
    try std.testing.expectEqual(73728, F[2].re);
    try std.testing.expectEqual(-9472, F[2].im);
    try std.testing.expectEqual(455168, F[3].re);
    try std.testing.expectEqual(-23040, F[3].im);
    try std.testing.expectEqual(1978496, F[4].re);
    try std.testing.expectEqual(-318336, F[4].im);
    try std.testing.expectEqual(23808, F[5].re);
    try std.testing.expectEqual(-7424, F[5].im);
    try std.testing.expectEqual(387456, F[6].re);
    try std.testing.expectEqual(33408, F[6].im);
    try std.testing.expectEqual(3743488, F[7].re);
    try std.testing.expectEqual(590592, F[7].im);
    try std.testing.expectEqual(33038464, F[8].re);
    try std.testing.expectEqual(12031360, F[8].im);
    try std.testing.expectEqual(78762752, F[9].re);
    try std.testing.expectEqual(-1556736, F[9].im);
    try std.testing.expectEqual(2912640, F[10].re);
    try std.testing.expectEqual(1066624, F[10].im);
    try std.testing.expectEqual(28104448, F[11].re);
    try std.testing.expectEqual(20824832, F[11].im);
    try std.testing.expectEqual(262698368, F[12].re);
    try std.testing.expectEqual(200420480, F[12].im);
    try std.testing.expectEqual(1236836096, F[13].re);
    try std.testing.expectEqual(1074371328, F[13].im);
    try std.testing.expectEqual(336448000, F[14].re);
    try std.testing.expectEqual(-2425088, F[14].im);
    try std.testing.expectEqual(68830464, F[15].re);
    try std.testing.expectEqual(46128384, F[15].im);
    try std.testing.expectEqual(706515200, F[16].re);
    try std.testing.expectEqual(551330560, F[16].im);
    try std.testing.expectEqual(3672819200, F[17].re);
    try std.testing.expectEqual(3366986240, F[17].im);
    try std.testing.expectEqual(7032608640, F[18].re);
    try std.testing.expectEqual(6507790720, F[18].im);
    try std.testing.expectEqual(942063360, F[19].re);
    try std.testing.expectEqual(43946240, F[19].im);

    trmm(.row_major, .right, .lower, .conj_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(512, F[0].re);
    try std.testing.expectEqual(-512, F[0].im);
    try std.testing.expectEqual(117760, F[1].re);
    try std.testing.expectEqual(-74752, F[1].im);
    try std.testing.expectEqual(4036096, F[2].re);
    try std.testing.expectEqual(-621056, F[2].im);
    try std.testing.expectEqual(40187904, F[3].re);
    try std.testing.expectEqual(-2615296, F[3].im);
    try std.testing.expectEqual(248699904, F[4].re);
    try std.testing.expectEqual(-35152896, F[4].im);
    try std.testing.expectEqual(95232, F[5].re);
    try std.testing.expectEqual(-29696, F[5].im);
    try std.testing.expectEqual(11420160, F[6].re);
    try std.testing.expectEqual(757248, F[6].im);
    try std.testing.expectEqual(214306816, F[7].re);
    try std.testing.expectEqual(31987712, F[7].im);
    try std.testing.expectEqual(2808325120, F[8].re);
    try std.testing.expectEqual(958702592, F[8].im);
    try std.testing.expectEqual(11428464640, F[9].re);
    try std.testing.expectEqual(1055987712, F[9].im);
    try std.testing.expectEqual(11650560, F[10].re);
    try std.testing.expectEqual(4266496, F[10].im);
    try std.testing.expectEqual(856827904, F[11].re);
    try std.testing.expectEqual(608694272, F[11].im);
    try std.testing.expectEqual(15137484800, F[12].re);
    try std.testing.expectEqual(11468388352, F[12].im);
    try std.testing.expectEqual(115011337216, F[13].re);
    try std.testing.expectEqual(97566848000, F[13].im);
    try std.testing.expectEqual(179267168256, F[14].re);
    try std.testing.expectEqual(123258004480, F[14].im);
    try std.testing.expectEqual(275321856, F[15].re);
    try std.testing.expectEqual(184513536, F[15].im);
    try std.testing.expectEqual(21434356736, F[16].re);
    try std.testing.expectEqual(16544336896, F[16].im);
    try std.testing.expectEqual(227927868416, F[17].re);
    try std.testing.expectEqual(203576800256, F[17].im);
    try std.testing.expectEqual(851369422336, F[18].re);
    try std.testing.expectEqual(777457798656, F[18].im);
    try std.testing.expectEqual(1175191228416, F[19].re);
    try std.testing.expectEqual(991297140736, F[19].im);

    trmm(.col_major, .right, .lower, .conj_trans, .non_unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(2048, F[0].re);
    try std.testing.expectEqual(-2048, F[0].im);
    try std.testing.expectEqual(471040, F[1].re);
    try std.testing.expectEqual(-299008, F[1].im);
    try std.testing.expectEqual(16144384, F[2].re);
    try std.testing.expectEqual(-2484224, F[2].im);
    try std.testing.expectEqual(160751616, F[3].re);
    try std.testing.expectEqual(-10461184, F[3].im);
    try std.testing.expectEqual(6963601408, F[4].re);
    try std.testing.expectEqual(-984285184, F[4].im);
    try std.testing.expectEqual(3608576, F[5].re);
    try std.testing.expectEqual(-1429504, F[5].im);
    try std.testing.expectEqual(352053248, F[6].re);
    try std.testing.expectEqual(16234496, F[6].im);
    try std.testing.expectEqual(6322094080, F[7].re);
    try std.testing.expectEqual(874733568, F[7].im);
    try std.testing.expectEqual(153991309312, F[8].re);
    try std.testing.expectEqual(48727635968, F[8].im);
    try std.testing.expectEqual(594284621824, F[9].re);
    try std.testing.expectEqual(54909513728, F[9].im);
    try std.testing.expectEqual(1019707392, F[10].re);
    try std.testing.expectEqual(238637056, F[10].im);
    try std.testing.expectEqual(51895123968, F[11].re);
    try std.testing.expectEqual(32644325376, F[11].im);
    try std.testing.expectEqual(1316668256256, F[12].re);
    try std.testing.expectEqual(924019347456, F[12].im);
    try std.testing.expectEqual(9380860960768, F[13].re);
    try std.testing.expectEqual(7474213494784, F[13].im);
    try std.testing.expectEqual(13625432922112, F[14].re);
    try std.testing.expectEqual(9367864588288, F[14].im);
    try std.testing.expectEqual(77264875520, F[15].re);
    try std.testing.expectEqual(49219620864, F[15].im);
    try std.testing.expectEqual(3532881971200, F[16].re);
    try std.testing.expectEqual(2628020787200, F[16].im);
    try std.testing.expectEqual(32679407861760, F[17].re);
    try std.testing.expectEqual(28226384445440, F[17].im);
    try std.testing.expectEqual(99479552256000, F[18].re);
    try std.testing.expectEqual(87606694082560, F[18].im);
    try std.testing.expectEqual(117601934295040, F[19].re);
    try std.testing.expectEqual(99182224015360, F[19].im);

    trmm(.row_major, .right, .lower, .conj_trans, .unit, m, n, beta, E.ptr, n, F.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(8192, F[0].re);
    try std.testing.expectEqual(0, F[0].im);
    try std.testing.expectEqual(1589248, F[1].re);
    try std.testing.expectEqual(294912, F[1].im);
    try std.testing.expectEqual(59957248, F[2].re);
    try std.testing.expectEqual(12877824, F[2].im);
    try std.testing.expectEqual(1536983040, F[3].re);
    try std.testing.expectEqual(101253120, F[3].im);
    try std.testing.expectEqual(32854835200, F[4].re);
    try std.testing.expectEqual(10699325440, F[4].im);
    try std.testing.expectEqual(10076160, F[5].re);
    try std.testing.expectEqual(4358144, F[5].im);
    try std.testing.expectEqual(758243328, F[6].re);
    try std.testing.expectEqual(702267392, F[6].im);
    try std.testing.expectEqual(27952054272, F[7].re);
    try std.testing.expectEqual(15110012928, F[7].im);
    try std.testing.expectEqual(689888690176, F[8].re);
    try std.testing.expectEqual(469431164928, F[8].im);
    try std.testing.expectEqual(16474832371712, F[9].re);
    try std.testing.expectEqual(6058025369600, F[9].im);
    try std.testing.expectEqual(1562140672, F[10].re);
    try std.testing.expectEqual(2516688896, F[10].im);
    try std.testing.expectEqual(62974574592, F[11].re);
    try std.testing.expectEqual(174806188032, F[11].im);
    try std.testing.expectEqual(3321130893312, F[12].re);
    try std.testing.expectEqual(6058802855936, F[12].im);
    try std.testing.expectEqual(102207539085312, F[13].re);
    try std.testing.expectEqual(102474628825088, F[13].im);
    try std.testing.expectEqual(1034863694807040, F[14].re);
    try std.testing.expectEqual(851413616631808, F[14].im);
    try std.testing.expectEqual(56090509312, F[15].re);
    try std.testing.expectEqual(252968992768, F[15].im);
    try std.testing.expectEqual(3664079380480, F[16].re);
    try std.testing.expectEqual(13503076417536, F[16].im);
    try std.testing.expectEqual(181884035973120, F[17].re);
    try std.testing.expectEqual(250122245718016, F[17].im);
    try std.testing.expectEqual(2621844008468480, F[18].re);
    try std.testing.expectEqual(2588327642013696, F[18].im);
    try std.testing.expectEqual(12910765823426560, F[19].re);
    try std.testing.expectEqual(11676038594953216, F[19].im);

    trmm(.col_major, .right, .lower, .conj_trans, .unit, m, n, beta, E.ptr, n, F.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(16384, F[0].re);
    try std.testing.expectEqual(16384, F[0].im);
    try std.testing.expectEqual(2588672, F[1].re);
    try std.testing.expectEqual(3768320, F[1].im);
    try std.testing.expectEqual(94158848, F[2].re);
    try std.testing.expectEqual(145670144, F[2].im);
    try std.testing.expectEqual(2871459840, F[3].re);
    try std.testing.expectEqual(3276472320, F[3].im);
    try std.testing.expectEqual(44311085056, F[4].re);
    try std.testing.expectEqual(87108321280, F[4].im);
    try std.testing.expectEqual(24150016, F[5].re);
    try std.testing.expectEqual(31227904, F[5].im);
    try std.testing.expectEqual(591609856, F[6].re);
    try std.testing.expectEqual(3024044032, F[6].im);
    try std.testing.expectEqual(37979947008, F[7].re);
    try std.testing.expectEqual(86934159360, F[7].im);
    try std.testing.expectEqual(1492269875200, F[8].re);
    try std.testing.expectEqual(2661018124288, F[8].im);
    try std.testing.expectEqual(20833955512320, F[9].re);
    try std.testing.expectEqual(45065858482176, F[9].im);
    try std.testing.expectEqual(23074177024, F[10].re);
    try std.testing.expectEqual(30784749568, F[10].im);
    try std.testing.expectEqual(689246306304, F[11].re);
    try std.testing.expectEqual(960296976384, F[11].im);
    try std.testing.expectEqual(34341196922880, F[12].re);
    try std.testing.expectEqual(45433188450304, F[12].im);
    try std.testing.expectEqual(922056821506048, F[13].re);
    try std.testing.expectEqual(748613918130176, F[13].im);
    try std.testing.expectEqual(367015892303872, F[14].re);
    try std.testing.expectEqual(3772721045127168, F[14].im);
    try std.testing.expectEqual(4163684892672, F[15].re);
    try std.testing.expectEqual(10952846049280, F[15].im);
    try std.testing.expectEqual(288719992373248, F[16].re);
    try std.testing.expectEqual(547632382984192, F[16].im);
    try std.testing.expectEqual(9028617084469248, F[17].re);
    try std.testing.expectEqual(9425464571789312, F[17].im);
    try std.testing.expectEqual(82856253574791170, F[18].re);
    try std.testing.expectEqual(78533611981094910, F[18].im);
    try std.testing.expectEqual(2478868993998848, F[19].re);
    try std.testing.expectEqual(49204941153042430, F[19].im);
}
