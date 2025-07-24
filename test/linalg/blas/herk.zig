const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const herk = zml.linalg.blas.herk;

test herk {
    const a = std.testing.allocator;

    const n = 5;
    const k = 3;
    const alpha: f64 = 2;
    const beta: f64 = 3;

    const A = try a.alloc(cf64, n * k);
    defer a.free(A);
    const B = try a.alloc(cf64, n * n);
    defer a.free(B);

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

    herk(.row_major, .upper, .no_trans, n, k, alpha, A.ptr, k, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(59, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(134, B[1].re);
    try std.testing.expectEqual(6, B[1].im);
    try std.testing.expectEqual(209, B[2].re);
    try std.testing.expectEqual(9, B[2].im);
    try std.testing.expectEqual(284, B[3].re);
    try std.testing.expectEqual(12, B[3].im);
    try std.testing.expectEqual(359, B[4].re);
    try std.testing.expectEqual(15, B[4].im);
    try std.testing.expectEqual(6, B[5].re);
    try std.testing.expectEqual(6, B[5].im);
    try std.testing.expectEqual(329, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(512, B[7].re);
    try std.testing.expectEqual(24, B[7].im);
    try std.testing.expectEqual(695, B[8].re);
    try std.testing.expectEqual(27, B[8].im);
    try std.testing.expectEqual(878, B[9].re);
    try std.testing.expectEqual(30, B[9].im);
    try std.testing.expectEqual(11, B[10].re);
    try std.testing.expectEqual(11, B[10].im);
    try std.testing.expectEqual(12, B[11].re);
    try std.testing.expectEqual(12, B[11].im);
    try std.testing.expectEqual(815, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(1106, B[13].re);
    try std.testing.expectEqual(42, B[13].im);
    try std.testing.expectEqual(1397, B[14].re);
    try std.testing.expectEqual(45, B[14].im);
    try std.testing.expectEqual(16, B[15].re);
    try std.testing.expectEqual(16, B[15].im);
    try std.testing.expectEqual(17, B[16].re);
    try std.testing.expectEqual(17, B[16].im);
    try std.testing.expectEqual(18, B[17].re);
    try std.testing.expectEqual(18, B[17].im);
    try std.testing.expectEqual(1517, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(1916, B[19].re);
    try std.testing.expectEqual(60, B[19].im);
    try std.testing.expectEqual(21, B[20].re);
    try std.testing.expectEqual(21, B[20].im);
    try std.testing.expectEqual(22, B[21].re);
    try std.testing.expectEqual(22, B[21].im);
    try std.testing.expectEqual(23, B[22].re);
    try std.testing.expectEqual(23, B[22].im);
    try std.testing.expectEqual(24, B[23].re);
    try std.testing.expectEqual(24, B[23].im);
    try std.testing.expectEqual(2435, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.col_major, .upper, .no_trans, n, k, alpha, A.ptr, n, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(809, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(134, B[1].re);
    try std.testing.expectEqual(6, B[1].im);
    try std.testing.expectEqual(209, B[2].re);
    try std.testing.expectEqual(9, B[2].im);
    try std.testing.expectEqual(284, B[3].re);
    try std.testing.expectEqual(12, B[3].im);
    try std.testing.expectEqual(359, B[4].re);
    try std.testing.expectEqual(15, B[4].im);
    try std.testing.expectEqual(722, B[5].re);
    try std.testing.expectEqual(18, B[5].im);
    try std.testing.expectEqual(1775, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(512, B[7].re);
    try std.testing.expectEqual(24, B[7].im);
    try std.testing.expectEqual(695, B[8].re);
    try std.testing.expectEqual(27, B[8].im);
    try std.testing.expectEqual(878, B[9].re);
    try std.testing.expectEqual(30, B[9].im);
    try std.testing.expectEqual(809, B[10].re);
    try std.testing.expectEqual(33, B[10].im);
    try std.testing.expectEqual(908, B[11].re);
    try std.testing.expectEqual(36, B[11].im);
    try std.testing.expectEqual(3413, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(1106, B[13].re);
    try std.testing.expectEqual(42, B[13].im);
    try std.testing.expectEqual(1397, B[14].re);
    try std.testing.expectEqual(45, B[14].im);
    try std.testing.expectEqual(896, B[15].re);
    try std.testing.expectEqual(48, B[15].im);
    try std.testing.expectEqual(1007, B[16].re);
    try std.testing.expectEqual(51, B[16].im);
    try std.testing.expectEqual(1118, B[17].re);
    try std.testing.expectEqual(54, B[17].im);
    try std.testing.expectEqual(5723, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(1916, B[19].re);
    try std.testing.expectEqual(60, B[19].im);
    try std.testing.expectEqual(983, B[20].re);
    try std.testing.expectEqual(63, B[20].im);
    try std.testing.expectEqual(1106, B[21].re);
    try std.testing.expectEqual(66, B[21].im);
    try std.testing.expectEqual(1229, B[22].re);
    try std.testing.expectEqual(69, B[22].im);
    try std.testing.expectEqual(1352, B[23].re);
    try std.testing.expectEqual(72, B[23].im);
    try std.testing.expectEqual(8705, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.row_major, .upper, .conj_trans, n, k, alpha, A.ptr, n, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(3059, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(1106, B[1].re);
    try std.testing.expectEqual(18, B[1].im);
    try std.testing.expectEqual(1403, B[2].re);
    try std.testing.expectEqual(27, B[2].im);
    try std.testing.expectEqual(1700, B[3].re);
    try std.testing.expectEqual(36, B[3].im);
    try std.testing.expectEqual(1997, B[4].re);
    try std.testing.expectEqual(45, B[4].im);
    try std.testing.expectEqual(722, B[5].re);
    try std.testing.expectEqual(18, B[5].im);
    try std.testing.expectEqual(6113, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(2408, B[7].re);
    try std.testing.expectEqual(72, B[7].im);
    try std.testing.expectEqual(3041, B[8].re);
    try std.testing.expectEqual(81, B[8].im);
    try std.testing.expectEqual(3674, B[9].re);
    try std.testing.expectEqual(90, B[9].im);
    try std.testing.expectEqual(809, B[10].re);
    try std.testing.expectEqual(33, B[10].im);
    try std.testing.expectEqual(908, B[11].re);
    try std.testing.expectEqual(36, B[11].im);
    try std.testing.expectEqual(11207, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(4382, B[13].re);
    try std.testing.expectEqual(126, B[13].im);
    try std.testing.expectEqual(5351, B[14].re);
    try std.testing.expectEqual(135, B[14].im);
    try std.testing.expectEqual(896, B[15].re);
    try std.testing.expectEqual(48, B[15].im);
    try std.testing.expectEqual(1007, B[16].re);
    try std.testing.expectEqual(51, B[16].im);
    try std.testing.expectEqual(1118, B[17].re);
    try std.testing.expectEqual(54, B[17].im);
    try std.testing.expectEqual(18341, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(7028, B[19].re);
    try std.testing.expectEqual(180, B[19].im);
    try std.testing.expectEqual(983, B[20].re);
    try std.testing.expectEqual(63, B[20].im);
    try std.testing.expectEqual(1106, B[21].re);
    try std.testing.expectEqual(66, B[21].im);
    try std.testing.expectEqual(1229, B[22].re);
    try std.testing.expectEqual(69, B[22].im);
    try std.testing.expectEqual(1352, B[23].re);
    try std.testing.expectEqual(72, B[23].im);
    try std.testing.expectEqual(27515, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.col_major, .upper, .conj_trans, n, k, alpha, A.ptr, k, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(9233, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(1106, B[1].re);
    try std.testing.expectEqual(18, B[1].im);
    try std.testing.expectEqual(1403, B[2].re);
    try std.testing.expectEqual(27, B[2].im);
    try std.testing.expectEqual(1700, B[3].re);
    try std.testing.expectEqual(36, B[3].im);
    try std.testing.expectEqual(1997, B[4].re);
    try std.testing.expectEqual(45, B[4].im);
    try std.testing.expectEqual(2294, B[5].re);
    try std.testing.expectEqual(54, B[5].im);
    try std.testing.expectEqual(18647, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(2408, B[7].re);
    try std.testing.expectEqual(72, B[7].im);
    try std.testing.expectEqual(3041, B[8].re);
    try std.testing.expectEqual(81, B[8].im);
    try std.testing.expectEqual(3674, B[9].re);
    try std.testing.expectEqual(90, B[9].im);
    try std.testing.expectEqual(2627, B[10].re);
    try std.testing.expectEqual(99, B[10].im);
    try std.testing.expectEqual(3212, B[11].re);
    try std.testing.expectEqual(108, B[11].im);
    try std.testing.expectEqual(34397, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(4382, B[13].re);
    try std.testing.expectEqual(126, B[13].im);
    try std.testing.expectEqual(5351, B[14].re);
    try std.testing.expectEqual(135, B[14].im);
    try std.testing.expectEqual(2960, B[15].re);
    try std.testing.expectEqual(144, B[15].im);
    try std.testing.expectEqual(3689, B[16].re);
    try std.testing.expectEqual(153, B[16].im);
    try std.testing.expectEqual(4418, B[17].re);
    try std.testing.expectEqual(162, B[17].im);
    try std.testing.expectEqual(56483, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(7028, B[19].re);
    try std.testing.expectEqual(180, B[19].im);
    try std.testing.expectEqual(3293, B[20].re);
    try std.testing.expectEqual(189, B[20].im);
    try std.testing.expectEqual(4166, B[21].re);
    try std.testing.expectEqual(198, B[21].im);
    try std.testing.expectEqual(5039, B[22].re);
    try std.testing.expectEqual(207, B[22].im);
    try std.testing.expectEqual(5912, B[23].re);
    try std.testing.expectEqual(216, B[23].im);
    try std.testing.expectEqual(84905, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.row_major, .lower, .no_trans, n, k, alpha, A.ptr, k, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(27755, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(1106, B[1].re);
    try std.testing.expectEqual(18, B[1].im);
    try std.testing.expectEqual(1403, B[2].re);
    try std.testing.expectEqual(27, B[2].im);
    try std.testing.expectEqual(1700, B[3].re);
    try std.testing.expectEqual(36, B[3].im);
    try std.testing.expectEqual(1997, B[4].re);
    try std.testing.expectEqual(45, B[4].im);
    try std.testing.expectEqual(7010, B[5].re);
    try std.testing.expectEqual(162, B[5].im);
    try std.testing.expectEqual(56249, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(2408, B[7].re);
    try std.testing.expectEqual(72, B[7].im);
    try std.testing.expectEqual(3041, B[8].re);
    try std.testing.expectEqual(81, B[8].im);
    try std.testing.expectEqual(3674, B[9].re);
    try std.testing.expectEqual(90, B[9].im);
    try std.testing.expectEqual(8081, B[10].re);
    try std.testing.expectEqual(297, B[10].im);
    try std.testing.expectEqual(10124, B[11].re);
    try std.testing.expectEqual(324, B[11].im);
    try std.testing.expectEqual(103967, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(4382, B[13].re);
    try std.testing.expectEqual(126, B[13].im);
    try std.testing.expectEqual(5351, B[14].re);
    try std.testing.expectEqual(135, B[14].im);
    try std.testing.expectEqual(9152, B[15].re);
    try std.testing.expectEqual(432, B[15].im);
    try std.testing.expectEqual(11735, B[16].re);
    try std.testing.expectEqual(459, B[16].im);
    try std.testing.expectEqual(14318, B[17].re);
    try std.testing.expectEqual(486, B[17].im);
    try std.testing.expectEqual(170909, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(7028, B[19].re);
    try std.testing.expectEqual(180, B[19].im);
    try std.testing.expectEqual(10223, B[20].re);
    try std.testing.expectEqual(567, B[20].im);
    try std.testing.expectEqual(13346, B[21].re);
    try std.testing.expectEqual(594, B[21].im);
    try std.testing.expectEqual(16469, B[22].re);
    try std.testing.expectEqual(621, B[22].im);
    try std.testing.expectEqual(19592, B[23].re);
    try std.testing.expectEqual(648, B[23].im);
    try std.testing.expectEqual(257075, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.col_major, .lower, .no_trans, n, k, alpha, A.ptr, n, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(83897, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(4022, B[1].re);
    try std.testing.expectEqual(54, B[1].im);
    try std.testing.expectEqual(4985, B[2].re);
    try std.testing.expectEqual(81, B[2].im);
    try std.testing.expectEqual(5948, B[3].re);
    try std.testing.expectEqual(108, B[3].im);
    try std.testing.expectEqual(6911, B[4].re);
    try std.testing.expectEqual(135, B[4].im);
    try std.testing.expectEqual(7010, B[5].re);
    try std.testing.expectEqual(162, B[5].im);
    try std.testing.expectEqual(169535, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(8096, B[7].re);
    try std.testing.expectEqual(216, B[7].im);
    try std.testing.expectEqual(10079, B[8].re);
    try std.testing.expectEqual(243, B[8].im);
    try std.testing.expectEqual(12062, B[9].re);
    try std.testing.expectEqual(270, B[9].im);
    try std.testing.expectEqual(8081, B[10].re);
    try std.testing.expectEqual(297, B[10].im);
    try std.testing.expectEqual(10124, B[11].re);
    try std.testing.expectEqual(324, B[11].im);
    try std.testing.expectEqual(312869, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(14210, B[13].re);
    try std.testing.expectEqual(378, B[13].im);
    try std.testing.expectEqual(17213, B[14].re);
    try std.testing.expectEqual(405, B[14].im);
    try std.testing.expectEqual(9152, B[15].re);
    try std.testing.expectEqual(432, B[15].im);
    try std.testing.expectEqual(11735, B[16].re);
    try std.testing.expectEqual(459, B[16].im);
    try std.testing.expectEqual(14318, B[17].re);
    try std.testing.expectEqual(486, B[17].im);
    try std.testing.expectEqual(513899, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(22364, B[19].re);
    try std.testing.expectEqual(540, B[19].im);
    try std.testing.expectEqual(10223, B[20].re);
    try std.testing.expectEqual(567, B[20].im);
    try std.testing.expectEqual(13346, B[21].re);
    try std.testing.expectEqual(594, B[21].im);
    try std.testing.expectEqual(16469, B[22].re);
    try std.testing.expectEqual(621, B[22].im);
    try std.testing.expectEqual(19592, B[23].re);
    try std.testing.expectEqual(648, B[23].im);
    try std.testing.expectEqual(772625, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.row_major, .lower, .conj_trans, n, k, alpha, A.ptr, n, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(252323, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(4022, B[1].re);
    try std.testing.expectEqual(54, B[1].im);
    try std.testing.expectEqual(4985, B[2].re);
    try std.testing.expectEqual(81, B[2].im);
    try std.testing.expectEqual(5948, B[3].re);
    try std.testing.expectEqual(108, B[3].im);
    try std.testing.expectEqual(6911, B[4].re);
    try std.testing.expectEqual(135, B[4].im);
    try std.testing.expectEqual(21734, B[5].re);
    try std.testing.expectEqual(486, B[5].im);
    try std.testing.expectEqual(509393, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(8096, B[7].re);
    try std.testing.expectEqual(216, B[7].im);
    try std.testing.expectEqual(10079, B[8].re);
    try std.testing.expectEqual(243, B[8].im);
    try std.testing.expectEqual(12062, B[9].re);
    try std.testing.expectEqual(270, B[9].im);
    try std.testing.expectEqual(25019, B[10].re);
    try std.testing.expectEqual(891, B[10].im);
    try std.testing.expectEqual(31244, B[11].re);
    try std.testing.expectEqual(972, B[11].im);
    try std.testing.expectEqual(939575, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(14210, B[13].re);
    try std.testing.expectEqual(378, B[13].im);
    try std.testing.expectEqual(17213, B[14].re);
    try std.testing.expectEqual(405, B[14].im);
    try std.testing.expectEqual(28304, B[15].re);
    try std.testing.expectEqual(1296, B[15].im);
    try std.testing.expectEqual(36161, B[16].re);
    try std.testing.expectEqual(1377, B[16].im);
    try std.testing.expectEqual(44018, B[17].re);
    try std.testing.expectEqual(1458, B[17].im);
    try std.testing.expectEqual(1542869, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(22364, B[19].re);
    try std.testing.expectEqual(540, B[19].im);
    try std.testing.expectEqual(31589, B[20].re);
    try std.testing.expectEqual(1701, B[20].im);
    try std.testing.expectEqual(41078, B[21].re);
    try std.testing.expectEqual(1782, B[21].im);
    try std.testing.expectEqual(50567, B[22].re);
    try std.testing.expectEqual(1863, B[22].im);
    try std.testing.expectEqual(60056, B[23].re);
    try std.testing.expectEqual(1944, B[23].im);
    try std.testing.expectEqual(2319275, B[24].re);
    try std.testing.expectEqual(0, B[24].im);

    herk(.col_major, .lower, .conj_trans, n, k, alpha, A.ptr, k, beta, B.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(757025, B[0].re);
    try std.testing.expectEqual(0, B[0].im);
    try std.testing.expectEqual(12194, B[1].re);
    try std.testing.expectEqual(162, B[1].im);
    try std.testing.expectEqual(15155, B[2].re);
    try std.testing.expectEqual(243, B[2].im);
    try std.testing.expectEqual(18116, B[3].re);
    try std.testing.expectEqual(324, B[3].im);
    try std.testing.expectEqual(21077, B[4].re);
    try std.testing.expectEqual(405, B[4].im);
    try std.testing.expectEqual(21734, B[5].re);
    try std.testing.expectEqual(486, B[5].im);
    try std.testing.expectEqual(1528487, B[6].re);
    try std.testing.expectEqual(0, B[6].im);
    try std.testing.expectEqual(24776, B[7].re);
    try std.testing.expectEqual(648, B[7].im);
    try std.testing.expectEqual(30905, B[8].re);
    try std.testing.expectEqual(729, B[8].im);
    try std.testing.expectEqual(37034, B[9].re);
    try std.testing.expectEqual(810, B[9].im);
    try std.testing.expectEqual(25019, B[10].re);
    try std.testing.expectEqual(891, B[10].im);
    try std.testing.expectEqual(31244, B[11].re);
    try std.testing.expectEqual(972, B[11].im);
    try std.testing.expectEqual(2819501, B[12].re);
    try std.testing.expectEqual(0, B[12].im);
    try std.testing.expectEqual(43694, B[13].re);
    try std.testing.expectEqual(1134, B[13].im);
    try std.testing.expectEqual(52991, B[14].re);
    try std.testing.expectEqual(1215, B[14].im);
    try std.testing.expectEqual(28304, B[15].re);
    try std.testing.expectEqual(1296, B[15].im);
    try std.testing.expectEqual(36161, B[16].re);
    try std.testing.expectEqual(1377, B[16].im);
    try std.testing.expectEqual(44018, B[17].re);
    try std.testing.expectEqual(1458, B[17].im);
    try std.testing.expectEqual(4630067, B[18].re);
    try std.testing.expectEqual(0, B[18].im);
    try std.testing.expectEqual(68948, B[19].re);
    try std.testing.expectEqual(1620, B[19].im);
    try std.testing.expectEqual(31589, B[20].re);
    try std.testing.expectEqual(1701, B[20].im);
    try std.testing.expectEqual(41078, B[21].re);
    try std.testing.expectEqual(1782, B[21].im);
    try std.testing.expectEqual(50567, B[22].re);
    try std.testing.expectEqual(1863, B[22].im);
    try std.testing.expectEqual(60056, B[23].re);
    try std.testing.expectEqual(1944, B[23].im);
    try std.testing.expectEqual(6960185, B[24].re);
    try std.testing.expectEqual(0, B[24].im);
}
