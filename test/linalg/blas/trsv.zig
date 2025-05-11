const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const trsv = zml.linalg.blas.trsv;

test trsv {
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

    trsv(f64, .RowMajor, .Upper, .NoTrans, .NonUnit, n, A.ptr, n, x1.ptr, 2);

    try std.testing.expectApproxEqRel(0, x1[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, x1[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, x1[4], 0.0000001);
    try std.testing.expectApproxEqRel(0, x1[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x1[8], 0.0000001);

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

    trsv(f64, .RowMajor, .Upper, .NoTrans, .Unit, n, A.ptr, n, x2.ptr, -2);

    try std.testing.expectApproxEqRel(15264, x2[8], 0.0000001);
    try std.testing.expectApproxEqRel(-9360, x2[6], 0.0000001);
    try std.testing.expectApproxEqRel(1272, x2[4], 0.0000001);
    try std.testing.expectApproxEqRel(-96, x2[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x2[0], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, A.ptr, n, x3.ptr, 2);

    try std.testing.expectApproxEqRel(-1.0364372469635623, x3[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044537, x3[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522267, x3[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789475, x3[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x3[8], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Upper, .NoTrans, .Unit, n, A.ptr, n, x4.ptr, -2);

    try std.testing.expectApproxEqRel(111104, x4[8], 0.0000001);
    try std.testing.expectApproxEqRel(-21848, x4[6], 0.0000001);
    try std.testing.expectApproxEqRel(1976, x4[4], 0.0000001);
    try std.testing.expectApproxEqRel(-116, x4[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x4[0], 0.0000001);

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

    trsv(f64, .RowMajor, .Upper, .Trans, .NonUnit, n, A.ptr, n, x5.ptr, 2);

    try std.testing.expectApproxEqRel(1, x5[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, x5[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, x5[4], 0.0000001);
    try std.testing.expectApproxEqRel(0, x5[6], 0.0000001);
    try std.testing.expectApproxEqRel(0, x5[8], 0.0000001);

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

    trsv(f64, .RowMajor, .Upper, .Trans, .Unit, n, A.ptr, n, x6.ptr, -2);

    try std.testing.expectApproxEqRel(1, x6[8], 0.0000001);
    try std.testing.expectApproxEqRel(0, x6[6], 0.0000001);
    try std.testing.expectApproxEqRel(0, x6[4], 0.0000001);
    try std.testing.expectApproxEqRel(0, x6[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, x6[0], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Upper, .Trans, .NonUnit, n, A.ptr, n, x7.ptr, 2);

    try std.testing.expectApproxEqRel(1, x7[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.5714285714285714, x7[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.08791208791208795, x7[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.037015615962984395, x7[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.020728744939271238, x7[8], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Upper, .Trans, .Unit, n, A.ptr, n, x8.ptr, -2);

    try std.testing.expectApproxEqRel(1, x8[8], 0.0000001);
    try std.testing.expectApproxEqRel(-4, x8[6], 0.0000001);
    try std.testing.expectApproxEqRel(40, x8[4], 0.0000001);
    try std.testing.expectApproxEqRel(-664, x8[2], 0.0000001);
    try std.testing.expectApproxEqRel(15088, x8[0], 0.0000001);

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

    trsv(f64, .RowMajor, .Lower, .NoTrans, .NonUnit, n, A.ptr, n, x9.ptr, 2);

    try std.testing.expectApproxEqRel(1, x9[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.5714285714285714, x9[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.08791208791208795, x9[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.037015615962984395, x9[6], 0.0000001);
    try std.testing.expectApproxEqRel(-0.020728744939271238, x9[8], 0.0000001);

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

    trsv(f64, .RowMajor, .Lower, .NoTrans, .Unit, n, A.ptr, n, x10.ptr, -2);

    try std.testing.expectApproxEqRel(1, x10[8], 0.0000001);
    try std.testing.expectApproxEqRel(-4, x10[6], 0.0000001);
    try std.testing.expectApproxEqRel(40, x10[4], 0.0000001);
    try std.testing.expectApproxEqRel(-664, x10[2], 0.0000001);
    try std.testing.expectApproxEqRel(15088, x10[0], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, A.ptr, n, x11.ptr, 2);

    try std.testing.expectApproxEqRel(1, x11[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, x11[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, x11[4], 0.0000001);
    try std.testing.expectApproxEqRel(0, x11[6], 0.0000001);
    try std.testing.expectApproxEqRel(0, x11[8], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Lower, .NoTrans, .Unit, n, A.ptr, n, x12.ptr, -2);

    try std.testing.expectApproxEqRel(1, x12[8], 0.0000001);
    try std.testing.expectApproxEqRel(0, x12[6], 0.0000001);
    try std.testing.expectApproxEqRel(0, x12[4], 0.0000001);
    try std.testing.expectApproxEqRel(0, x12[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, x12[0], 0.0000001);

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

    trsv(f64, .RowMajor, .Lower, .Trans, .NonUnit, n, A.ptr, n, x13.ptr, 2);

    try std.testing.expectApproxEqRel(-1.0364372469635623, x13[0], 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044537, x13[2], 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522267, x13[4], 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789475, x13[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x13[8], 0.0000001);

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

    trsv(f64, .RowMajor, .Lower, .Trans, .Unit, n, A.ptr, n, x14.ptr, -2);

    try std.testing.expectApproxEqRel(111104, x14[8], 0.0000001);
    try std.testing.expectApproxEqRel(-21848, x14[6], 0.0000001);
    try std.testing.expectApproxEqRel(1976, x14[4], 0.0000001);
    try std.testing.expectApproxEqRel(-116, x14[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x14[0], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Lower, .Trans, .NonUnit, n, A.ptr, n, x15.ptr, 2);

    try std.testing.expectApproxEqRel(0, x15[0], 0.0000001);
    try std.testing.expectApproxEqRel(0, x15[2], 0.0000001);
    try std.testing.expectApproxEqRel(0, x15[4], 0.0000001);
    try std.testing.expectApproxEqRel(0, x15[6], 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x15[8], 0.0000001);

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

    trsv(f64, .ColumnMajor, .Lower, .Trans, .Unit, n, A.ptr, n, x16.ptr, -2);

    try std.testing.expectApproxEqRel(15264, x16[8], 0.0000001);
    try std.testing.expectApproxEqRel(-9360, x16[6], 0.0000001);
    try std.testing.expectApproxEqRel(1272, x16[4], 0.0000001);
    try std.testing.expectApproxEqRel(-96, x16[2], 0.0000001);
    try std.testing.expectApproxEqRel(5, x16[0], 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .NoTrans, .NonUnit, n, B.ptr, n, x17.ptr, 2);

    try std.testing.expectApproxEqRel(0, x17[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x17[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x17[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .NoTrans, .Unit, n, B.ptr, n, x18.ptr, -2);

    try std.testing.expectApproxEqRel(-59241, x18[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-87681, x18[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(40906, x18[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(3678, x18[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2797, x18[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(2541, x18[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x18[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-196, x18[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x18[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x18[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .NoTrans, .NonUnit, n, B.ptr, n, x19.ptr, 2);

    try std.testing.expectApproxEqRel(-1.0364372469635625, x19[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044532, x19[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522266, x19[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789478, x19[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x19[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x19[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .NoTrans, .Unit, n, B.ptr, n, x20.ptr, -2);

    try std.testing.expectApproxEqRel(-434745, x20[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-611985, x20[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(95114, x20[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(8142, x20[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4317, x20[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(3949, x20[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x20[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-236, x20[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x20[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x20[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .ConjNoTrans, .NonUnit, n, B.ptr, n, x21.ptr, 2);

    try std.testing.expectApproxEqRel(0, x21[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x21[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x21[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .ConjNoTrans, .Unit, n, B.ptr, n, x22.ptr, -2);

    try std.testing.expectApproxEqRel(-87681, x22[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-59241, x22[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(3678, x22[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(40906, x22[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(2541, x22[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2797, x22[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-196, x22[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x22[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x22[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x22[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .ConjNoTrans, .NonUnit, n, B.ptr, n, x23.ptr, 2);

    try std.testing.expectApproxEqRel(0, x23[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.0364372469635628, x23[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044532, x23[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522267, x23[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789477, x23[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x23[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x23[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .ConjNoTrans, .Unit, n, B.ptr, n, x24.ptr, -2);

    try std.testing.expectApproxEqRel(-611985, x24[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-434745, x24[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(8142, x24[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(95114, x24[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(3949, x24[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4317, x24[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-236, x24[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x24[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x24[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x24[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .Trans, .NonUnit, n, B.ptr, n, x25.ptr, 2);

    try std.testing.expectApproxEqRel(1, x25[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x25[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .Trans, .Unit, n, B.ptr, n, x26.ptr, -2);

    try std.testing.expectApproxEqRel(1, x26[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x26[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x26[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x26[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-29, x26[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3, x26[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(332, x26[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(444, x26[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(2595, x26[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15045, x26[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .Trans, .NonUnit, n, B.ptr, n, x27.ptr, 2);

    try std.testing.expectApproxEqRel(1, x27[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5714285714285714, x27[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08791208791208795, x27[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03701561596298439, x27[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.020728744939271238, x27[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x27[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .Trans, .Unit, n, B.ptr, n, x28.ptr, -2);

    try std.testing.expectApproxEqRel(1, x28[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x28[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x28[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10, x28[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-141, x28[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(77, x28[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(3724, x28[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1260, x28[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-54381, x28[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-118005, x28[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .ConjTrans, .NonUnit, n, B.ptr, n, x29.ptr, 2);

    try std.testing.expectApproxEqRel(0, x29[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x29[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x29[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Upper, .ConjTrans, .Unit, n, B.ptr, n, x30.ptr, -2);

    try std.testing.expectApproxEqRel(1, x30[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x30[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x30[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x30[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3, x30[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29, x30[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(444, x30[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(332, x30[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-15045, x30[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2595, x30[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .ConjTrans, .NonUnit, n, B.ptr, n, x31.ptr, 2);

    try std.testing.expectApproxEqRel(0, x31[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x31[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5714285714285714, x31[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08791208791208795, x31[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03701561596298439, x31[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x31[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.020728744939271238, x31[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Upper, .ConjTrans, .Unit, n, B.ptr, n, x32.ptr, -2);

    try std.testing.expectApproxEqRel(1, x32[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x32[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10, x32[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x32[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(77, x32[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-141, x32[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1260, x32[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(3724, x32[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-118005, x32[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-54381, x32[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .NoTrans, .NonUnit, n, B.ptr, n, x33.ptr, 2);

    try std.testing.expectApproxEqRel(1, x33[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5714285714285714, x33[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08791208791208795, x33[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03701561596298439, x33[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.020728744939271238, x33[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x33[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .NoTrans, .Unit, n, B.ptr, n, x34.ptr, -2);

    try std.testing.expectApproxEqRel(1, x34[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x34[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x34[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-10, x34[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-141, x34[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(77, x34[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(3724, x34[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(1260, x34[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-54381, x34[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-118005, x34[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .NoTrans, .NonUnit, n, B.ptr, n, x35.ptr, 2);

    try std.testing.expectApproxEqRel(1, x35[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x35[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .NoTrans, .Unit, n, B.ptr, n, x36.ptr, -2);

    try std.testing.expectApproxEqRel(1, x36[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x36[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(2, x36[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x36[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-29, x36[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-3, x36[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(332, x36[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(444, x36[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(2595, x36[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-15045, x36[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .ConjNoTrans, .NonUnit, n, B.ptr, n, x37.ptr, 2);

    try std.testing.expectApproxEqRel(0, x37[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x37[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.5714285714285714, x37[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.08791208791208795, x37[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.03701561596298439, x37[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x37[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.020728744939271238, x37[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .ConjNoTrans, .Unit, n, B.ptr, n, x38.ptr, -2);

    try std.testing.expectApproxEqRel(1, x38[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x38[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-10, x38[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x38[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(77, x38[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-141, x38[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(1260, x38[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(3724, x38[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-118005, x38[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-54381, x38[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .ConjNoTrans, .NonUnit, n, B.ptr, n, x39.ptr, 2);

    try std.testing.expectApproxEqRel(0, x39[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x39[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x39[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .ConjNoTrans, .Unit, n, B.ptr, n, x40.ptr, -2);

    try std.testing.expectApproxEqRel(1, x40[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(1, x40[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2, x40[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(2, x40[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-3, x40[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-29, x40[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(444, x40[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(332, x40[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-15045, x40[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(2595, x40[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .Trans, .NonUnit, n, B.ptr, n, x41.ptr, 2);

    try std.testing.expectApproxEqRel(-1.0364372469635628, x41[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044532, x41[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522267, x41[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789477, x41[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x41[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x41[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .Trans, .Unit, n, B.ptr, n, x42.ptr, -2);

    try std.testing.expectApproxEqRel(-434745, x42[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-611985, x42[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(95114, x42[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(8142, x42[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-4317, x42[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(3949, x42[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x42[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-236, x42[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x42[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x42[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .Trans, .NonUnit, n, B.ptr, n, x43.ptr, 2);

    try std.testing.expectApproxEqRel(0, x43[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x43[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x43[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .Trans, .Unit, n, B.ptr, n, x44.ptr, -2);

    try std.testing.expectApproxEqRel(-59241, x44[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-87681, x44[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(40906, x44[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(3678, x44[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(-2797, x44[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(2541, x44[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(4, x44[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-196, x44[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x44[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x44[0].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .ConjTrans, .NonUnit, n, B.ptr, n, x45.ptr, 2);

    try std.testing.expectApproxEqRel(0, x45[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(-1.0364372469635628, x45[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.12955465587044532, x45[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.06477732793522267, x45[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(-0.04210526315789477, x45[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x45[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x45[8].im, 0.0000001);

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

    trsv(cf64, .RowMajor, .Lower, .ConjTrans, .Unit, n, B.ptr, n, x46.ptr, -2);

    try std.testing.expectApproxEqRel(-611985, x46[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-434745, x46[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(8142, x46[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(95114, x46[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(3949, x46[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-4317, x46[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-236, x46[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x46[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x46[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x46[0].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .ConjTrans, .NonUnit, n, B.ptr, n, x47.ptr, 2);

    try std.testing.expectApproxEqRel(0, x47[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[0].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(0, x47[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(0.2, x47[8].im, 0.0000001);

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

    trsv(cf64, .ColumnMajor, .Lower, .ConjTrans, .Unit, n, B.ptr, n, x48.ptr, -2);

    try std.testing.expectApproxEqRel(-87681, x48[8].re, 0.0000001);
    try std.testing.expectApproxEqRel(-59241, x48[8].im, 0.0000001);
    try std.testing.expectApproxEqRel(3678, x48[6].re, 0.0000001);
    try std.testing.expectApproxEqRel(40906, x48[6].im, 0.0000001);
    try std.testing.expectApproxEqRel(2541, x48[4].re, 0.0000001);
    try std.testing.expectApproxEqRel(-2797, x48[4].im, 0.0000001);
    try std.testing.expectApproxEqRel(-196, x48[2].re, 0.0000001);
    try std.testing.expectApproxEqRel(4, x48[2].im, 0.0000001);
    try std.testing.expectApproxEqRel(5, x48[0].re, 0.0000001);
    try std.testing.expectApproxEqRel(5, x48[0].im, 0.0000001);
}
