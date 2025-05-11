const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const tpmv = zml.linalg.blas.tpmv;

test tpmv {
    const a = std.testing.allocator;

    const n = 5;

    const A = try a.alloc(f64, n * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);

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
        21,
        22,
        23,
        24,
        25,
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

    tpmv(f64, .RowMajor, .Upper, .NoTrans, .NonUnit, n, A.ptr, x1.ptr, 2);

    try std.testing.expectEqual(55, x1[0]);
    try std.testing.expectEqual(110, x1[2]);
    try std.testing.expectEqual(134, x1[4]);
    try std.testing.expectEqual(122, x1[6]);
    try std.testing.expectEqual(75, x1[8]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);

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

    tpmv(f64, .RowMajor, .Upper, .NoTrans, .Unit, n, A.ptr, x2.ptr, -2);

    try std.testing.expectEqual(55, x2[8]);
    try std.testing.expectEqual(100, x2[6]);
    try std.testing.expectEqual(107, x2[4]);
    try std.testing.expectEqual(74, x2[2]);
    try std.testing.expectEqual(5, x2[0]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);

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

    tpmv(f64, .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, A.ptr, x3.ptr, 2);

    try std.testing.expectEqual(100, x3[0]);
    try std.testing.expectEqual(113, x3[2]);
    try std.testing.expectEqual(119, x3[4]);
    try std.testing.expectEqual(110, x3[6]);
    try std.testing.expectEqual(75, x3[8]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);

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

    tpmv(f64, .ColumnMajor, .Upper, .NoTrans, .Unit, n, A.ptr, x4.ptr, -2);

    try std.testing.expectEqual(100, x4[8]);
    try std.testing.expectEqual(109, x4[6]);
    try std.testing.expectEqual(104, x4[4]);
    try std.testing.expectEqual(74, x4[2]);
    try std.testing.expectEqual(5, x4[0]);

    const x5 = try a.alloc(f64, 2 * n);
    defer a.free(x5);

    @memcpy(x5.ptr, &[_]f64{
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

    tpmv(f64, .RowMajor, .Upper, .Trans, .NonUnit, n, A.ptr, x5.ptr, 2);

    try std.testing.expectEqual(1, x5[0]);
    try std.testing.expectEqual(14, x5[2]);
    try std.testing.expectEqual(47, x5[4]);
    try std.testing.expectEqual(105, x5[6]);
    try std.testing.expectEqual(190, x5[8]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);

    @memcpy(x6.ptr, &[_]f64{
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

    tpmv(f64, .RowMajor, .Upper, .Trans, .Unit, n, A.ptr, x6.ptr, -2);

    try std.testing.expectEqual(1, x6[8]);
    try std.testing.expectEqual(4, x6[6]);
    try std.testing.expectEqual(20, x6[4]);
    try std.testing.expectEqual(57, x6[2]);
    try std.testing.expectEqual(120, x6[0]);

    const x7 = try a.alloc(f64, 2 * n);
    defer a.free(x7);

    @memcpy(x7.ptr, &[_]f64{
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

    tpmv(f64, .ColumnMajor, .Upper, .Trans, .NonUnit, n, A.ptr, x7.ptr, 2);

    try std.testing.expectEqual(1, x7[0]);
    try std.testing.expectEqual(8, x7[2]);
    try std.testing.expectEqual(32, x7[4]);
    try std.testing.expectEqual(90, x7[6]);
    try std.testing.expectEqual(205, x7[8]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);

    @memcpy(x8.ptr, &[_]f64{
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

    tpmv(f64, .ColumnMajor, .Upper, .Trans, .Unit, n, A.ptr, x8.ptr, -2);

    try std.testing.expectEqual(1, x8[8]);
    try std.testing.expectEqual(4, x8[6]);
    try std.testing.expectEqual(17, x8[4]);
    try std.testing.expectEqual(54, x8[2]);
    try std.testing.expectEqual(135, x8[0]);

    const x9 = try a.alloc(f64, 2 * n);
    defer a.free(x9);

    @memcpy(x9.ptr, &[_]f64{
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

    tpmv(f64, .RowMajor, .Lower, .NoTrans, .NonUnit, n, A.ptr, x9.ptr, 2);

    try std.testing.expectEqual(1, x9[0]);
    try std.testing.expectEqual(8, x9[2]);
    try std.testing.expectEqual(32, x9[4]);
    try std.testing.expectEqual(90, x9[6]);
    try std.testing.expectEqual(205, x9[8]);

    const x10 = try a.alloc(f64, 2 * n);
    defer a.free(x10);

    @memcpy(x10.ptr, &[_]f64{
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

    tpmv(f64, .RowMajor, .Lower, .NoTrans, .Unit, n, A.ptr, x10.ptr, -2);

    try std.testing.expectEqual(1, x10[8]);
    try std.testing.expectEqual(4, x10[6]);
    try std.testing.expectEqual(17, x10[4]);
    try std.testing.expectEqual(54, x10[2]);
    try std.testing.expectEqual(135, x10[0]);

    const x11 = try a.alloc(f64, 2 * n);
    defer a.free(x11);

    @memcpy(x11.ptr, &[_]f64{
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

    tpmv(f64, .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, A.ptr, x11.ptr, 2);

    try std.testing.expectEqual(1, x11[0]);
    try std.testing.expectEqual(14, x11[2]);
    try std.testing.expectEqual(47, x11[4]);
    try std.testing.expectEqual(105, x11[6]);
    try std.testing.expectEqual(190, x11[8]);

    const x12 = try a.alloc(f64, 2 * n);
    defer a.free(x12);

    @memcpy(x12.ptr, &[_]f64{
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

    tpmv(f64, .ColumnMajor, .Lower, .NoTrans, .Unit, n, A.ptr, x12.ptr, -2);

    try std.testing.expectEqual(1, x12[8]);
    try std.testing.expectEqual(4, x12[6]);
    try std.testing.expectEqual(20, x12[4]);
    try std.testing.expectEqual(57, x12[2]);
    try std.testing.expectEqual(120, x12[0]);

    const x13 = try a.alloc(f64, 2 * n);
    defer a.free(x13);

    @memcpy(x13.ptr, &[_]f64{
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

    tpmv(f64, .RowMajor, .Lower, .Trans, .NonUnit, n, A.ptr, x13.ptr, 2);

    try std.testing.expectEqual(100, x13[0]);
    try std.testing.expectEqual(113, x13[2]);
    try std.testing.expectEqual(119, x13[4]);
    try std.testing.expectEqual(110, x13[6]);
    try std.testing.expectEqual(75, x13[8]);

    const x14 = try a.alloc(f64, 2 * n);
    defer a.free(x14);

    @memcpy(x14.ptr, &[_]f64{
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

    tpmv(f64, .RowMajor, .Lower, .Trans, .Unit, n, A.ptr, x14.ptr, -2);

    try std.testing.expectEqual(100, x14[8]);
    try std.testing.expectEqual(109, x14[6]);
    try std.testing.expectEqual(104, x14[4]);
    try std.testing.expectEqual(74, x14[2]);
    try std.testing.expectEqual(5, x14[0]);

    const x15 = try a.alloc(f64, 2 * n);
    defer a.free(x15);

    @memcpy(x15.ptr, &[_]f64{
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

    tpmv(f64, .ColumnMajor, .Lower, .Trans, .NonUnit, n, A.ptr, x15.ptr, 2);

    try std.testing.expectEqual(55, x15[0]);
    try std.testing.expectEqual(110, x15[2]);
    try std.testing.expectEqual(134, x15[4]);
    try std.testing.expectEqual(122, x15[6]);
    try std.testing.expectEqual(75, x15[8]);

    const x16 = try a.alloc(f64, 2 * n);
    defer a.free(x16);

    @memcpy(x16.ptr, &[_]f64{
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

    tpmv(f64, .ColumnMajor, .Lower, .Trans, .Unit, n, A.ptr, x16.ptr, -2);

    try std.testing.expectEqual(55, x16[8]);
    try std.testing.expectEqual(100, x16[6]);
    try std.testing.expectEqual(107, x16[4]);
    try std.testing.expectEqual(74, x16[2]);
    try std.testing.expectEqual(5, x16[0]);

    const B = try a.alloc(cf64, n * n);
    defer a.free(B);
    const x17 = try a.alloc(cf64, 2 * n);
    defer a.free(x17);

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
    @memcpy(x17.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .NoTrans, .NonUnit, n, B.ptr, x17.ptr, 2);

    try std.testing.expectEqual(0, x17[0].re);
    try std.testing.expectEqual(110, x17[0].im);
    try std.testing.expectEqual(0, x17[2].re);
    try std.testing.expectEqual(220, x17[2].im);
    try std.testing.expectEqual(0, x17[4].re);
    try std.testing.expectEqual(268, x17[4].im);
    try std.testing.expectEqual(0, x17[6].re);
    try std.testing.expectEqual(244, x17[6].im);
    try std.testing.expectEqual(0, x17[8].re);
    try std.testing.expectEqual(150, x17[8].im);

    const x18 = try a.alloc(cf64, 2 * n);
    defer a.free(x18);

    @memcpy(x18.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .NoTrans, .Unit, n, B.ptr, x18.ptr, -2);

    try std.testing.expectEqual(1, x18[8].re);
    try std.testing.expectEqual(109, x18[8].im);
    try std.testing.expectEqual(2, x18[6].re);
    try std.testing.expectEqual(198, x18[6].im);
    try std.testing.expectEqual(3, x18[4].re);
    try std.testing.expectEqual(211, x18[4].im);
    try std.testing.expectEqual(4, x18[2].re);
    try std.testing.expectEqual(144, x18[2].im);
    try std.testing.expectEqual(5, x18[0].re);
    try std.testing.expectEqual(5, x18[0].im);

    const x19 = try a.alloc(cf64, 2 * n);
    defer a.free(x19);

    @memcpy(x19.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, B.ptr, x19.ptr, 2);

    try std.testing.expectEqual(0, x19[0].re);
    try std.testing.expectEqual(200, x19[0].im);
    try std.testing.expectEqual(0, x19[2].re);
    try std.testing.expectEqual(226, x19[2].im);
    try std.testing.expectEqual(0, x19[4].re);
    try std.testing.expectEqual(238, x19[4].im);
    try std.testing.expectEqual(0, x19[6].re);
    try std.testing.expectEqual(220, x19[6].im);
    try std.testing.expectEqual(0, x19[8].re);
    try std.testing.expectEqual(150, x19[8].im);

    const x20 = try a.alloc(cf64, 2 * n);
    defer a.free(x20);

    @memcpy(x20.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .NoTrans, .Unit, n, B.ptr, x20.ptr, -2);

    try std.testing.expectEqual(1, x20[8].re);
    try std.testing.expectEqual(199, x20[8].im);
    try std.testing.expectEqual(2, x20[6].re);
    try std.testing.expectEqual(216, x20[6].im);
    try std.testing.expectEqual(3, x20[4].re);
    try std.testing.expectEqual(205, x20[4].im);
    try std.testing.expectEqual(4, x20[2].re);
    try std.testing.expectEqual(144, x20[2].im);
    try std.testing.expectEqual(5, x20[0].re);
    try std.testing.expectEqual(5, x20[0].im);

    const x21 = try a.alloc(cf64, 2 * n);
    defer a.free(x21);

    @memcpy(x21.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .ConjNoTrans, .NonUnit, n, B.ptr, x21.ptr, 2);

    try std.testing.expectEqual(110, x21[0].re);
    try std.testing.expectEqual(0, x21[0].im);
    try std.testing.expectEqual(220, x21[2].re);
    try std.testing.expectEqual(0, x21[2].im);
    try std.testing.expectEqual(268, x21[4].re);
    try std.testing.expectEqual(0, x21[4].im);
    try std.testing.expectEqual(244, x21[6].re);
    try std.testing.expectEqual(0, x21[6].im);
    try std.testing.expectEqual(150, x21[8].re);
    try std.testing.expectEqual(0, x21[8].im);

    const x22 = try a.alloc(cf64, 2 * n);
    defer a.free(x22);

    @memcpy(x22.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .ConjNoTrans, .Unit, n, B.ptr, x22.ptr, -2);

    try std.testing.expectEqual(109, x22[8].re);
    try std.testing.expectEqual(1, x22[8].im);
    try std.testing.expectEqual(198, x22[6].re);
    try std.testing.expectEqual(2, x22[6].im);
    try std.testing.expectEqual(211, x22[4].re);
    try std.testing.expectEqual(3, x22[4].im);
    try std.testing.expectEqual(144, x22[2].re);
    try std.testing.expectEqual(4, x22[2].im);
    try std.testing.expectEqual(5, x22[0].re);
    try std.testing.expectEqual(5, x22[0].im);

    const x23 = try a.alloc(cf64, 2 * n);
    defer a.free(x23);

    @memcpy(x23.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .ConjNoTrans, .NonUnit, n, B.ptr, x23.ptr, 2);

    try std.testing.expectEqual(200, x23[0].re);
    try std.testing.expectEqual(0, x23[0].im);
    try std.testing.expectEqual(226, x23[2].re);
    try std.testing.expectEqual(0, x23[2].im);
    try std.testing.expectEqual(238, x23[4].re);
    try std.testing.expectEqual(0, x23[4].im);
    try std.testing.expectEqual(220, x23[6].re);
    try std.testing.expectEqual(0, x23[6].im);
    try std.testing.expectEqual(150, x23[8].re);
    try std.testing.expectEqual(0, x23[8].im);

    const x24 = try a.alloc(cf64, 2 * n);
    defer a.free(x24);

    @memcpy(x24.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .ConjNoTrans, .Unit, n, B.ptr, x24.ptr, -2);

    try std.testing.expectEqual(199, x24[8].re);
    try std.testing.expectEqual(1, x24[8].im);
    try std.testing.expectEqual(216, x24[6].re);
    try std.testing.expectEqual(2, x24[6].im);
    try std.testing.expectEqual(205, x24[4].re);
    try std.testing.expectEqual(3, x24[4].im);
    try std.testing.expectEqual(144, x24[2].re);
    try std.testing.expectEqual(4, x24[2].im);
    try std.testing.expectEqual(5, x24[0].re);
    try std.testing.expectEqual(5, x24[0].im);

    const x25 = try a.alloc(cf64, 2 * n);
    defer a.free(x25);

    @memcpy(x25.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .Trans, .NonUnit, n, B.ptr, x25.ptr, 2);

    try std.testing.expectEqual(0, x25[0].re);
    try std.testing.expectEqual(2, x25[0].im);
    try std.testing.expectEqual(0, x25[2].re);
    try std.testing.expectEqual(28, x25[2].im);
    try std.testing.expectEqual(0, x25[4].re);
    try std.testing.expectEqual(94, x25[4].im);
    try std.testing.expectEqual(0, x25[6].re);
    try std.testing.expectEqual(210, x25[6].im);
    try std.testing.expectEqual(0, x25[8].re);
    try std.testing.expectEqual(380, x25[8].im);

    const x26 = try a.alloc(cf64, 2 * n);
    defer a.free(x26);

    @memcpy(x26.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .Trans, .Unit, n, B.ptr, x26.ptr, -2);

    try std.testing.expectEqual(1, x26[8].re);
    try std.testing.expectEqual(1, x26[8].im);
    try std.testing.expectEqual(2, x26[6].re);
    try std.testing.expectEqual(6, x26[6].im);
    try std.testing.expectEqual(3, x26[4].re);
    try std.testing.expectEqual(37, x26[4].im);
    try std.testing.expectEqual(4, x26[2].re);
    try std.testing.expectEqual(110, x26[2].im);
    try std.testing.expectEqual(5, x26[0].re);
    try std.testing.expectEqual(235, x26[0].im);

    const x27 = try a.alloc(cf64, 2 * n);
    defer a.free(x27);

    @memcpy(x27.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .Trans, .NonUnit, n, B.ptr, x27.ptr, 2);

    try std.testing.expectEqual(0, x27[0].re);
    try std.testing.expectEqual(2, x27[0].im);
    try std.testing.expectEqual(0, x27[2].re);
    try std.testing.expectEqual(16, x27[2].im);
    try std.testing.expectEqual(0, x27[4].re);
    try std.testing.expectEqual(64, x27[4].im);
    try std.testing.expectEqual(0, x27[6].re);
    try std.testing.expectEqual(180, x27[6].im);
    try std.testing.expectEqual(0, x27[8].re);
    try std.testing.expectEqual(410, x27[8].im);

    const x28 = try a.alloc(cf64, 2 * n);
    defer a.free(x28);

    @memcpy(x28.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .Trans, .Unit, n, B.ptr, x28.ptr, -2);

    try std.testing.expectEqual(1, x28[8].re);
    try std.testing.expectEqual(1, x28[8].im);
    try std.testing.expectEqual(2, x28[6].re);
    try std.testing.expectEqual(6, x28[6].im);
    try std.testing.expectEqual(3, x28[4].re);
    try std.testing.expectEqual(31, x28[4].im);
    try std.testing.expectEqual(4, x28[2].re);
    try std.testing.expectEqual(104, x28[2].im);
    try std.testing.expectEqual(5, x28[0].re);
    try std.testing.expectEqual(265, x28[0].im);

    const x29 = try a.alloc(cf64, 2 * n);
    defer a.free(x29);

    @memcpy(x29.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .ConjTrans, .NonUnit, n, B.ptr, x29.ptr, 2);

    try std.testing.expectEqual(2, x29[0].re);
    try std.testing.expectEqual(0, x29[0].im);
    try std.testing.expectEqual(28, x29[2].re);
    try std.testing.expectEqual(0, x29[2].im);
    try std.testing.expectEqual(94, x29[4].re);
    try std.testing.expectEqual(0, x29[4].im);
    try std.testing.expectEqual(210, x29[6].re);
    try std.testing.expectEqual(0, x29[6].im);
    try std.testing.expectEqual(380, x29[8].re);
    try std.testing.expectEqual(0, x29[8].im);

    const x30 = try a.alloc(cf64, 2 * n);
    defer a.free(x30);

    @memcpy(x30.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Upper, .ConjTrans, .Unit, n, B.ptr, x30.ptr, -2);

    try std.testing.expectEqual(1, x30[8].re);
    try std.testing.expectEqual(1, x30[8].im);
    try std.testing.expectEqual(6, x30[6].re);
    try std.testing.expectEqual(2, x30[6].im);
    try std.testing.expectEqual(37, x30[4].re);
    try std.testing.expectEqual(3, x30[4].im);
    try std.testing.expectEqual(110, x30[2].re);
    try std.testing.expectEqual(4, x30[2].im);
    try std.testing.expectEqual(235, x30[0].re);
    try std.testing.expectEqual(5, x30[0].im);

    const x31 = try a.alloc(cf64, 2 * n);
    defer a.free(x31);

    @memcpy(x31.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .ConjTrans, .NonUnit, n, B.ptr, x31.ptr, 2);

    try std.testing.expectEqual(2, x31[0].re);
    try std.testing.expectEqual(0, x31[0].im);
    try std.testing.expectEqual(16, x31[2].re);
    try std.testing.expectEqual(0, x31[2].im);
    try std.testing.expectEqual(64, x31[4].re);
    try std.testing.expectEqual(0, x31[4].im);
    try std.testing.expectEqual(180, x31[6].re);
    try std.testing.expectEqual(0, x31[6].im);
    try std.testing.expectEqual(410, x31[8].re);
    try std.testing.expectEqual(0, x31[8].im);

    const x32 = try a.alloc(cf64, 2 * n);
    defer a.free(x32);

    @memcpy(x32.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Upper, .ConjTrans, .Unit, n, B.ptr, x32.ptr, -2);

    try std.testing.expectEqual(1, x32[8].re);
    try std.testing.expectEqual(1, x32[8].im);
    try std.testing.expectEqual(6, x32[6].re);
    try std.testing.expectEqual(2, x32[6].im);
    try std.testing.expectEqual(31, x32[4].re);
    try std.testing.expectEqual(3, x32[4].im);
    try std.testing.expectEqual(104, x32[2].re);
    try std.testing.expectEqual(4, x32[2].im);
    try std.testing.expectEqual(265, x32[0].re);
    try std.testing.expectEqual(5, x32[0].im);

    const x33 = try a.alloc(cf64, 2 * n);
    defer a.free(x33);

    @memcpy(x33.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .NoTrans, .NonUnit, n, B.ptr, x33.ptr, 2);

    try std.testing.expectEqual(0, x33[0].re);
    try std.testing.expectEqual(2, x33[0].im);
    try std.testing.expectEqual(0, x33[2].re);
    try std.testing.expectEqual(16, x33[2].im);
    try std.testing.expectEqual(0, x33[4].re);
    try std.testing.expectEqual(64, x33[4].im);
    try std.testing.expectEqual(0, x33[6].re);
    try std.testing.expectEqual(180, x33[6].im);
    try std.testing.expectEqual(0, x33[8].re);
    try std.testing.expectEqual(410, x33[8].im);

    const x34 = try a.alloc(cf64, 2 * n);
    defer a.free(x34);

    @memcpy(x34.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .NoTrans, .Unit, n, B.ptr, x34.ptr, -2);

    try std.testing.expectEqual(1, x34[8].re);
    try std.testing.expectEqual(1, x34[8].im);
    try std.testing.expectEqual(2, x34[6].re);
    try std.testing.expectEqual(6, x34[6].im);
    try std.testing.expectEqual(3, x34[4].re);
    try std.testing.expectEqual(31, x34[4].im);
    try std.testing.expectEqual(4, x34[2].re);
    try std.testing.expectEqual(104, x34[2].im);
    try std.testing.expectEqual(5, x34[0].re);
    try std.testing.expectEqual(265, x34[0].im);

    const x35 = try a.alloc(cf64, 2 * n);
    defer a.free(x35);

    @memcpy(x35.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, B.ptr, x35.ptr, 2);

    try std.testing.expectEqual(0, x35[0].re);
    try std.testing.expectEqual(2, x35[0].im);
    try std.testing.expectEqual(0, x35[2].re);
    try std.testing.expectEqual(28, x35[2].im);
    try std.testing.expectEqual(0, x35[4].re);
    try std.testing.expectEqual(94, x35[4].im);
    try std.testing.expectEqual(0, x35[6].re);
    try std.testing.expectEqual(210, x35[6].im);
    try std.testing.expectEqual(0, x35[8].re);
    try std.testing.expectEqual(380, x35[8].im);

    const x36 = try a.alloc(cf64, 2 * n);
    defer a.free(x36);

    @memcpy(x36.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .NoTrans, .Unit, n, B.ptr, x36.ptr, -2);

    try std.testing.expectEqual(1, x36[8].re);
    try std.testing.expectEqual(1, x36[8].im);
    try std.testing.expectEqual(2, x36[6].re);
    try std.testing.expectEqual(6, x36[6].im);
    try std.testing.expectEqual(3, x36[4].re);
    try std.testing.expectEqual(37, x36[4].im);
    try std.testing.expectEqual(4, x36[2].re);
    try std.testing.expectEqual(110, x36[2].im);
    try std.testing.expectEqual(5, x36[0].re);
    try std.testing.expectEqual(235, x36[0].im);

    const x37 = try a.alloc(cf64, 2 * n);
    defer a.free(x37);

    @memcpy(x37.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .ConjNoTrans, .NonUnit, n, B.ptr, x37.ptr, 2);

    try std.testing.expectEqual(2, x37[0].re);
    try std.testing.expectEqual(0, x37[0].im);
    try std.testing.expectEqual(16, x37[2].re);
    try std.testing.expectEqual(0, x37[2].im);
    try std.testing.expectEqual(64, x37[4].re);
    try std.testing.expectEqual(0, x37[4].im);
    try std.testing.expectEqual(180, x37[6].re);
    try std.testing.expectEqual(0, x37[6].im);
    try std.testing.expectEqual(410, x37[8].re);
    try std.testing.expectEqual(0, x37[8].im);

    const x38 = try a.alloc(cf64, 2 * n);
    defer a.free(x38);

    @memcpy(x38.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .ConjNoTrans, .Unit, n, B.ptr, x38.ptr, -2);

    try std.testing.expectEqual(1, x38[8].re);
    try std.testing.expectEqual(1, x38[8].im);
    try std.testing.expectEqual(6, x38[6].re);
    try std.testing.expectEqual(2, x38[6].im);
    try std.testing.expectEqual(31, x38[4].re);
    try std.testing.expectEqual(3, x38[4].im);
    try std.testing.expectEqual(104, x38[2].re);
    try std.testing.expectEqual(4, x38[2].im);
    try std.testing.expectEqual(265, x38[0].re);
    try std.testing.expectEqual(5, x38[0].im);

    const x39 = try a.alloc(cf64, 2 * n);
    defer a.free(x39);

    @memcpy(x39.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .ConjNoTrans, .NonUnit, n, B.ptr, x39.ptr, 2);

    try std.testing.expectEqual(2, x39[0].re);
    try std.testing.expectEqual(0, x39[0].im);
    try std.testing.expectEqual(28, x39[2].re);
    try std.testing.expectEqual(0, x39[2].im);
    try std.testing.expectEqual(94, x39[4].re);
    try std.testing.expectEqual(0, x39[4].im);
    try std.testing.expectEqual(210, x39[6].re);
    try std.testing.expectEqual(0, x39[6].im);
    try std.testing.expectEqual(380, x39[8].re);
    try std.testing.expectEqual(0, x39[8].im);

    const x40 = try a.alloc(cf64, 2 * n);
    defer a.free(x40);

    @memcpy(x40.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .ConjNoTrans, .Unit, n, B.ptr, x40.ptr, -2);

    try std.testing.expectEqual(1, x40[8].re);
    try std.testing.expectEqual(1, x40[8].im);
    try std.testing.expectEqual(6, x40[6].re);
    try std.testing.expectEqual(2, x40[6].im);
    try std.testing.expectEqual(37, x40[4].re);
    try std.testing.expectEqual(3, x40[4].im);
    try std.testing.expectEqual(110, x40[2].re);
    try std.testing.expectEqual(4, x40[2].im);
    try std.testing.expectEqual(235, x40[0].re);
    try std.testing.expectEqual(5, x40[0].im);

    const x41 = try a.alloc(cf64, 2 * n);
    defer a.free(x41);

    @memcpy(x41.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .Trans, .NonUnit, n, B.ptr, x41.ptr, 2);

    try std.testing.expectEqual(0, x41[0].re);
    try std.testing.expectEqual(200, x41[0].im);
    try std.testing.expectEqual(0, x41[2].re);
    try std.testing.expectEqual(226, x41[2].im);
    try std.testing.expectEqual(0, x41[4].re);
    try std.testing.expectEqual(238, x41[4].im);
    try std.testing.expectEqual(0, x41[6].re);
    try std.testing.expectEqual(220, x41[6].im);
    try std.testing.expectEqual(0, x41[8].re);
    try std.testing.expectEqual(150, x41[8].im);

    const x42 = try a.alloc(cf64, 2 * n);
    defer a.free(x42);

    @memcpy(x42.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .Trans, .Unit, n, B.ptr, x42.ptr, -2);

    try std.testing.expectEqual(1, x42[8].re);
    try std.testing.expectEqual(199, x42[8].im);
    try std.testing.expectEqual(2, x42[6].re);
    try std.testing.expectEqual(216, x42[6].im);
    try std.testing.expectEqual(3, x42[4].re);
    try std.testing.expectEqual(205, x42[4].im);
    try std.testing.expectEqual(4, x42[2].re);
    try std.testing.expectEqual(144, x42[2].im);
    try std.testing.expectEqual(5, x42[0].re);
    try std.testing.expectEqual(5, x42[0].im);

    const x43 = try a.alloc(cf64, 2 * n);
    defer a.free(x43);

    @memcpy(x43.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .Trans, .NonUnit, n, B.ptr, x43.ptr, 2);

    try std.testing.expectEqual(0, x43[0].re);
    try std.testing.expectEqual(110, x43[0].im);
    try std.testing.expectEqual(0, x43[2].re);
    try std.testing.expectEqual(220, x43[2].im);
    try std.testing.expectEqual(0, x43[4].re);
    try std.testing.expectEqual(268, x43[4].im);
    try std.testing.expectEqual(0, x43[6].re);
    try std.testing.expectEqual(244, x43[6].im);
    try std.testing.expectEqual(0, x43[8].re);
    try std.testing.expectEqual(150, x43[8].im);

    const x44 = try a.alloc(cf64, 2 * n);
    defer a.free(x44);

    @memcpy(x44.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .Trans, .Unit, n, B.ptr, x44.ptr, -2);

    try std.testing.expectEqual(1, x44[8].re);
    try std.testing.expectEqual(109, x44[8].im);
    try std.testing.expectEqual(2, x44[6].re);
    try std.testing.expectEqual(198, x44[6].im);
    try std.testing.expectEqual(3, x44[4].re);
    try std.testing.expectEqual(211, x44[4].im);
    try std.testing.expectEqual(4, x44[2].re);
    try std.testing.expectEqual(144, x44[2].im);
    try std.testing.expectEqual(5, x44[0].re);
    try std.testing.expectEqual(5, x44[0].im);

    const x45 = try a.alloc(cf64, 2 * n);
    defer a.free(x45);

    @memcpy(x45.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .ConjTrans, .NonUnit, n, B.ptr, x45.ptr, 2);

    try std.testing.expectEqual(200, x45[0].re);
    try std.testing.expectEqual(0, x45[0].im);
    try std.testing.expectEqual(226, x45[2].re);
    try std.testing.expectEqual(0, x45[2].im);
    try std.testing.expectEqual(238, x45[4].re);
    try std.testing.expectEqual(0, x45[4].im);
    try std.testing.expectEqual(220, x45[6].re);
    try std.testing.expectEqual(0, x45[6].im);
    try std.testing.expectEqual(150, x45[8].re);
    try std.testing.expectEqual(0, x45[8].im);

    const x46 = try a.alloc(cf64, 2 * n);
    defer a.free(x46);

    @memcpy(x46.ptr, &[_]cf64{
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

    tpmv(cf64, .RowMajor, .Lower, .ConjTrans, .Unit, n, B.ptr, x46.ptr, -2);

    try std.testing.expectEqual(199, x46[8].re);
    try std.testing.expectEqual(1, x46[8].im);
    try std.testing.expectEqual(216, x46[6].re);
    try std.testing.expectEqual(2, x46[6].im);
    try std.testing.expectEqual(205, x46[4].re);
    try std.testing.expectEqual(3, x46[4].im);
    try std.testing.expectEqual(144, x46[2].re);
    try std.testing.expectEqual(4, x46[2].im);
    try std.testing.expectEqual(5, x46[0].re);
    try std.testing.expectEqual(5, x46[0].im);

    const x47 = try a.alloc(cf64, 2 * n);
    defer a.free(x47);

    @memcpy(x47.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .ConjTrans, .NonUnit, n, B.ptr, x47.ptr, 2);

    try std.testing.expectEqual(110, x47[0].re);
    try std.testing.expectEqual(0, x47[0].im);
    try std.testing.expectEqual(220, x47[2].re);
    try std.testing.expectEqual(0, x47[2].im);
    try std.testing.expectEqual(268, x47[4].re);
    try std.testing.expectEqual(0, x47[4].im);
    try std.testing.expectEqual(244, x47[6].re);
    try std.testing.expectEqual(0, x47[6].im);
    try std.testing.expectEqual(150, x47[8].re);
    try std.testing.expectEqual(0, x47[8].im);

    const x48 = try a.alloc(cf64, 2 * n);
    defer a.free(x48);

    @memcpy(x48.ptr, &[_]cf64{
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

    tpmv(cf64, .ColumnMajor, .Lower, .ConjTrans, .Unit, n, B.ptr, x48.ptr, -2);

    try std.testing.expectEqual(109, x48[8].re);
    try std.testing.expectEqual(1, x48[8].im);
    try std.testing.expectEqual(198, x48[6].re);
    try std.testing.expectEqual(2, x48[6].im);
    try std.testing.expectEqual(211, x48[4].re);
    try std.testing.expectEqual(3, x48[4].im);
    try std.testing.expectEqual(144, x48[2].re);
    try std.testing.expectEqual(4, x48[2].im);
    try std.testing.expectEqual(5, x48[0].re);
    try std.testing.expectEqual(5, x48[0].im);
}
