const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const gbmv = zml.linalg.blas.gbmv;

test gbmv {
    const a = std.testing.allocator;

    const m = 4;
    const n = 5;
    const kl = 1;
    const ku = 2;
    const alpha: f64 = 2;
    const beta: f64 = 3;

    const A = try a.alloc(f64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, 2 * m);
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
    });

    gbmv(.row_major, .no_trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x1.ptr, 2, beta, y1.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(43, y1[0]);
    try std.testing.expectEqual(146, y1[2]);
    try std.testing.expectEqual(313, y1[4]);
    try std.testing.expectEqual(352, y1[6]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * m);
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
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    gbmv(.row_major, .no_trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x2.ptr, -2, beta, y2.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(43, y2[6]);
    try std.testing.expectEqual(146, y2[4]);
    try std.testing.expectEqual(313, y2[2]);
    try std.testing.expectEqual(352, y2[0]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * m);
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
    });

    gbmv(.col_major, .no_trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x3.ptr, 2, beta, y3.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(87, y3[0]);
    try std.testing.expectEqual(206, y3[2]);
    try std.testing.expectEqual(389, y3[4]);
    try std.testing.expectEqual(384, y3[6]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(f64, 2 * m);
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
        4,
        0,
        3,
        0,
        2,
        0,
        1,
        0,
    });

    gbmv(.col_major, .no_trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x4.ptr, -2, beta, y4.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(87, y4[6]);
    try std.testing.expectEqual(206, y4[4]);
    try std.testing.expectEqual(389, y4[2]);
    try std.testing.expectEqual(384, y4[0]);

    const x5 = try a.alloc(f64, 2 * m);
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

    gbmv(.row_major, .trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x5.ptr, 2, beta, y5.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(27, y5[0]);
    try std.testing.expectEqual(90, y5[2]);
    try std.testing.expectEqual(209, y5[4]);
    try std.testing.expectEqual(222, y5[6]);
    try std.testing.expectEqual(207, y5[8]);

    const x6 = try a.alloc(f64, 2 * m);
    defer a.free(x6);
    const y6 = try a.alloc(f64, 2 * n);
    defer a.free(y6);

    @memcpy(x6.ptr, &[_]f64{
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

    gbmv(.row_major, .trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x6.ptr, -2, beta, y6.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(27, y6[8]);
    try std.testing.expectEqual(90, y6[6]);
    try std.testing.expectEqual(209, y6[4]);
    try std.testing.expectEqual(222, y6[2]);
    try std.testing.expectEqual(207, y6[0]);

    const x7 = try a.alloc(f64, 2 * m);
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

    gbmv(.col_major, .trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x7.ptr, 2, beta, y7.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(25, y7[0]);
    try std.testing.expectEqual(94, y7[2]);
    try std.testing.expectEqual(229, y7[4]);
    try std.testing.expectEqual(268, y7[6]);
    try std.testing.expectEqual(261, y7[8]);

    const x8 = try a.alloc(f64, 2 * m);
    defer a.free(x8);
    const y8 = try a.alloc(f64, 2 * n);
    defer a.free(y8);

    @memcpy(x8.ptr, &[_]f64{
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

    gbmv(.col_major, .trans, m, n, kl, ku, alpha, A.ptr, kl + ku + 1, x8.ptr, -2, beta, y8.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(25, y8[8]);
    try std.testing.expectEqual(94, y8[6]);
    try std.testing.expectEqual(229, y8[4]);
    try std.testing.expectEqual(268, y8[2]);
    try std.testing.expectEqual(261, y8[0]);

    const gamma = cf64.init(2, 2);
    const delta = cf64.init(3, 3);

    const B = try a.alloc(cf64, m * n);
    defer a.free(B);
    const x9 = try a.alloc(cf64, 2 * n);
    defer a.free(x9);
    const y9 = try a.alloc(cf64, 2 * m);
    defer a.free(y9);

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
    });
    @memcpy(x9.ptr, &[_]cf64{
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
    @memcpy(y9.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
    });

    gbmv(.row_major, .no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x9.ptr, 2, delta, y9.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(-160, y9[0].re);
    try std.testing.expectEqual(130, y9[0].im);
    try std.testing.expectEqual(-560, y9[2].re);
    try std.testing.expectEqual(468, y9[2].im);
    try std.testing.expectEqual(-1216, y9[4].re);
    try std.testing.expectEqual(1066, y9[4].im);
    try std.testing.expectEqual(-1360, y9[6].re);
    try std.testing.expectEqual(1216, y9[6].im);

    const x10 = try a.alloc(cf64, 2 * n);
    defer a.free(x10);
    const y10 = try a.alloc(cf64, 2 * m);
    defer a.free(y10);

    @memcpy(x10.ptr, &[_]cf64{
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
    @memcpy(y10.ptr, &[_]cf64{
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    gbmv(.row_major, .no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x10.ptr, -2, delta, y10.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(-160, y10[6].re);
    try std.testing.expectEqual(130, y10[6].im);
    try std.testing.expectEqual(-560, y10[4].re);
    try std.testing.expectEqual(468, y10[4].im);
    try std.testing.expectEqual(-1216, y10[2].re);
    try std.testing.expectEqual(1066, y10[2].im);
    try std.testing.expectEqual(-1360, y10[0].re);
    try std.testing.expectEqual(1216, y10[0].im);

    const x11 = try a.alloc(cf64, 2 * n);
    defer a.free(x11);
    const y11 = try a.alloc(cf64, 2 * m);
    defer a.free(y11);

    @memcpy(x11.ptr, &[_]cf64{
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
    @memcpy(y11.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
    });

    gbmv(.col_major, .no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x11.ptr, 2, delta, y11.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(-336, y11[0].re);
    try std.testing.expectEqual(270, y11[0].im);
    try std.testing.expectEqual(-800, y11[2].re);
    try std.testing.expectEqual(676, y11[2].im);
    try std.testing.expectEqual(-1520, y11[4].re);
    try std.testing.expectEqual(1338, y11[4].im);
    try std.testing.expectEqual(-1488, y11[6].re);
    try std.testing.expectEqual(1332, y11[6].im);

    const x12 = try a.alloc(cf64, 2 * n);
    defer a.free(x12);
    const y12 = try a.alloc(cf64, 2 * m);
    defer a.free(y12);

    @memcpy(x12.ptr, &[_]cf64{
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
    @memcpy(y12.ptr, &[_]cf64{
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    gbmv(.col_major, .no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x12.ptr, -2, delta, y12.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(-336, y12[6].re);
    try std.testing.expectEqual(270, y12[6].im);
    try std.testing.expectEqual(-800, y12[4].re);
    try std.testing.expectEqual(676, y12[4].im);
    try std.testing.expectEqual(-1520, y12[2].re);
    try std.testing.expectEqual(1338, y12[2].im);
    try std.testing.expectEqual(-1488, y12[0].re);
    try std.testing.expectEqual(1332, y12[0].im);

    const x13 = try a.alloc(cf64, 2 * n);
    defer a.free(x13);
    const y13 = try a.alloc(cf64, 2 * m);
    defer a.free(y13);

    @memcpy(x13.ptr, &[_]cf64{
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
    @memcpy(y13.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
    });

    gbmv(.row_major, .conj_no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x13.ptr, 2, delta, y13.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(124, y13[0].re);
    try std.testing.expectEqual(166, y13[0].im);
    try std.testing.expectEqual(456, y13[2].re);
    try std.testing.expectEqual(572, y13[2].im);
    try std.testing.expectEqual(1048, y13[4].re);
    try std.testing.expectEqual(1234, y13[4].im);
    try std.testing.expectEqual(1192, y13[6].re);
    try std.testing.expectEqual(1384, y13[6].im);

    const x14 = try a.alloc(cf64, 2 * n);
    defer a.free(x14);
    const y14 = try a.alloc(cf64, 2 * m);
    defer a.free(y14);

    @memcpy(x14.ptr, &[_]cf64{
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
    @memcpy(y14.ptr, &[_]cf64{
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    gbmv(.row_major, .conj_no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x14.ptr, -2, delta, y14.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(124, y14[6].re);
    try std.testing.expectEqual(166, y14[6].im);
    try std.testing.expectEqual(456, y14[4].re);
    try std.testing.expectEqual(572, y14[4].im);
    try std.testing.expectEqual(1048, y14[2].re);
    try std.testing.expectEqual(1234, y14[2].im);
    try std.testing.expectEqual(1192, y14[0].re);
    try std.testing.expectEqual(1384, y14[0].im);

    const x15 = try a.alloc(cf64, 2 * n);
    defer a.free(x15);
    const y15 = try a.alloc(cf64, 2 * m);
    defer a.free(y15);

    @memcpy(x15.ptr, &[_]cf64{
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
    @memcpy(y15.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
    });

    gbmv(.col_major, .conj_no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x15.ptr, 2, delta, y15.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(264, y15[0].re);
    try std.testing.expectEqual(342, y15[0].im);
    try std.testing.expectEqual(664, y15[2].re);
    try std.testing.expectEqual(812, y15[2].im);
    try std.testing.expectEqual(1320, y15[4].re);
    try std.testing.expectEqual(1538, y15[4].im);
    try std.testing.expectEqual(1308, y15[6].re);
    try std.testing.expectEqual(1512, y15[6].im);

    const x16 = try a.alloc(cf64, 2 * n);
    defer a.free(x16);
    const y16 = try a.alloc(cf64, 2 * m);
    defer a.free(y16);

    @memcpy(x16.ptr, &[_]cf64{
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
    @memcpy(y16.ptr, &[_]cf64{
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    gbmv(.col_major, .conj_no_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x16.ptr, -2, delta, y16.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(264, y16[6].re);
    try std.testing.expectEqual(342, y16[6].im);
    try std.testing.expectEqual(664, y16[4].re);
    try std.testing.expectEqual(812, y16[4].im);
    try std.testing.expectEqual(1320, y16[2].re);
    try std.testing.expectEqual(1538, y16[2].im);
    try std.testing.expectEqual(1308, y16[0].re);
    try std.testing.expectEqual(1512, y16[0].im);

    const x17 = try a.alloc(cf64, 2 * m);
    defer a.free(x17);
    const y17 = try a.alloc(cf64, 2 * n);
    defer a.free(y17);

    @memcpy(x17.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
    });
    @memcpy(y17.ptr, &[_]cf64{
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

    gbmv(.row_major, .trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x17.ptr, 2, delta, y17.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(-96, y17[0].re);
    try std.testing.expectEqual(74, y17[0].im);
    try std.testing.expectEqual(-336, y17[2].re);
    try std.testing.expectEqual(276, y17[2].im);
    try std.testing.expectEqual(-800, y17[4].re);
    try std.testing.expectEqual(682, y17[4].im);
    try std.testing.expectEqual(-840, y17[6].re);
    try std.testing.expectEqual(732, y17[6].im);
    try std.testing.expectEqual(-768, y17[8].re);
    try std.testing.expectEqual(690, y17[8].im);

    const x18 = try a.alloc(cf64, 2 * m);
    defer a.free(x18);
    const y18 = try a.alloc(cf64, 2 * n);
    defer a.free(y18);

    @memcpy(x18.ptr, &[_]cf64{
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y18.ptr, &[_]cf64{
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

    gbmv(.row_major, .trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x18.ptr, -2, delta, y18.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(-96, y18[8].re);
    try std.testing.expectEqual(74, y18[8].im);
    try std.testing.expectEqual(-336, y18[6].re);
    try std.testing.expectEqual(276, y18[6].im);
    try std.testing.expectEqual(-800, y18[4].re);
    try std.testing.expectEqual(682, y18[4].im);
    try std.testing.expectEqual(-840, y18[2].re);
    try std.testing.expectEqual(732, y18[2].im);
    try std.testing.expectEqual(-768, y18[0].re);
    try std.testing.expectEqual(690, y18[0].im);

    const x19 = try a.alloc(cf64, 2 * m);
    defer a.free(x19);
    const y19 = try a.alloc(cf64, 2 * n);
    defer a.free(y19);

    @memcpy(x19.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
    });
    @memcpy(y19.ptr, &[_]cf64{
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

    gbmv(.col_major, .trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x19.ptr, 2, delta, y19.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(-88, y19[0].re);
    try std.testing.expectEqual(66, y19[0].im);
    try std.testing.expectEqual(-352, y19[2].re);
    try std.testing.expectEqual(280, y19[2].im);
    try std.testing.expectEqual(-880, y19[4].re);
    try std.testing.expectEqual(730, y19[4].im);
    try std.testing.expectEqual(-1024, y19[6].re);
    try std.testing.expectEqual(880, y19[6].im);
    try std.testing.expectEqual(-984, y19[8].re);
    try std.testing.expectEqual(874, y19[8].im);

    const x20 = try a.alloc(cf64, 2 * m);
    defer a.free(x20);
    const y20 = try a.alloc(cf64, 2 * n);
    defer a.free(y20);

    @memcpy(x20.ptr, &[_]cf64{
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y20.ptr, &[_]cf64{
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

    gbmv(.col_major, .trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x20.ptr, -2, delta, y20.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(-88, y20[8].re);
    try std.testing.expectEqual(66, y20[8].im);
    try std.testing.expectEqual(-352, y20[6].re);
    try std.testing.expectEqual(280, y20[6].im);
    try std.testing.expectEqual(-880, y20[4].re);
    try std.testing.expectEqual(730, y20[4].im);
    try std.testing.expectEqual(-1024, y20[2].re);
    try std.testing.expectEqual(880, y20[2].im);
    try std.testing.expectEqual(-984, y20[0].re);
    try std.testing.expectEqual(874, y20[0].im);

    const x21 = try a.alloc(cf64, 2 * m);
    defer a.free(x21);
    const y21 = try a.alloc(cf64, 2 * n);
    defer a.free(y21);

    @memcpy(x21.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
    });
    @memcpy(y21.ptr, &[_]cf64{
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

    gbmv(.row_major, .conj_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x21.ptr, 2, delta, y21.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(68, y21[0].re);
    try std.testing.expectEqual(102, y21[0].im);
    try std.testing.expectEqual(264, y21[2].re);
    try std.testing.expectEqual(348, y21[2].im);
    try std.testing.expectEqual(664, y21[4].re);
    try std.testing.expectEqual(818, y21[4].im);
    try std.testing.expectEqual(708, y21[6].re);
    try std.testing.expectEqual(864, y21[6].im);
    try std.testing.expectEqual(660, y21[8].re);
    try std.testing.expectEqual(798, y21[8].im);

    const x22 = try a.alloc(cf64, 2 * m);
    defer a.free(x22);
    const y22 = try a.alloc(cf64, 2 * n);
    defer a.free(y22);

    @memcpy(x22.ptr, &[_]cf64{
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y22.ptr, &[_]cf64{
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

    gbmv(.row_major, .conj_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x22.ptr, -2, delta, y22.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(68, y22[8].re);
    try std.testing.expectEqual(102, y22[8].im);
    try std.testing.expectEqual(264, y22[6].re);
    try std.testing.expectEqual(348, y22[6].im);
    try std.testing.expectEqual(664, y22[4].re);
    try std.testing.expectEqual(818, y22[4].im);
    try std.testing.expectEqual(708, y22[2].re);
    try std.testing.expectEqual(864, y22[2].im);
    try std.testing.expectEqual(660, y22[0].re);
    try std.testing.expectEqual(798, y22[0].im);

    const x23 = try a.alloc(cf64, 2 * m);
    defer a.free(x23);
    const y23 = try a.alloc(cf64, 2 * n);
    defer a.free(y23);

    @memcpy(x23.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
    });
    @memcpy(y23.ptr, &[_]cf64{
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

    gbmv(.col_major, .conj_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x23.ptr, 2, delta, y23.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(60, y23[0].re);
    try std.testing.expectEqual(94, y23[0].im);
    try std.testing.expectEqual(268, y23[2].re);
    try std.testing.expectEqual(364, y23[2].im);
    try std.testing.expectEqual(712, y23[4].re);
    try std.testing.expectEqual(898, y23[4].im);
    try std.testing.expectEqual(856, y23[6].re);
    try std.testing.expectEqual(1048, y23[6].im);
    try std.testing.expectEqual(844, y23[8].re);
    try std.testing.expectEqual(1014, y23[8].im);

    const x24 = try a.alloc(cf64, 2 * m);
    defer a.free(x24);
    const y24 = try a.alloc(cf64, 2 * n);
    defer a.free(y24);

    @memcpy(x24.ptr, &[_]cf64{
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y24.ptr, &[_]cf64{
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

    gbmv(.col_major, .conj_trans, m, n, kl, ku, gamma, B.ptr, kl + ku + 1, x24.ptr, -2, delta, y24.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(60, y24[8].re);
    try std.testing.expectEqual(94, y24[8].im);
    try std.testing.expectEqual(268, y24[6].re);
    try std.testing.expectEqual(364, y24[6].im);
    try std.testing.expectEqual(712, y24[4].re);
    try std.testing.expectEqual(898, y24[4].im);
    try std.testing.expectEqual(856, y24[2].re);
    try std.testing.expectEqual(1048, y24[2].im);
    try std.testing.expectEqual(844, y24[0].re);
    try std.testing.expectEqual(1014, y24[0].im);
}
