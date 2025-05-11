const std = @import("std");
const zml = @import("zml");
const sbmv = zml.linalg.blas.sbmv;

test sbmv {
    const a = std.testing.allocator;

    const n = 5;
    const k = 1;
    const alpha = 2;
    const beta = 3;

    const A = try a.alloc(f64, (1 + 2 * k) * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * n);
    defer a.free(x1);
    const y1 = try a.alloc(f64, 2 * n);
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
        5,
        0,
    });

    sbmv(f64, .RowMajor, .Upper, n, k, alpha, A.ptr, k + 1, x1.ptr, 2, beta, y1.ptr, 2);

    try std.testing.expectEqual(13, y1[0]);
    try std.testing.expectEqual(46, y1[2]);
    try std.testing.expectEqual(103, y1[4]);
    try std.testing.expectEqual(184, y1[6]);
    try std.testing.expectEqual(169, y1[8]);

    const x2 = try a.alloc(f64, 2 * n);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * n);
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

    sbmv(f64, .RowMajor, .Upper, n, k, alpha, A.ptr, k + 1, x2.ptr, -2, beta, y2.ptr, -2);

    try std.testing.expectEqual(13, y2[8]);
    try std.testing.expectEqual(46, y2[6]);
    try std.testing.expectEqual(103, y2[4]);
    try std.testing.expectEqual(184, y2[2]);
    try std.testing.expectEqual(169, y2[0]);

    const x3 = try a.alloc(f64, 2 * n);
    defer a.free(x3);
    const y3 = try a.alloc(f64, 2 * n);
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
        5,
        0,
    });

    sbmv(f64, .ColumnMajor, .Upper, n, k, alpha, A.ptr, k + 1, x3.ptr, 2, beta, y3.ptr, 2);

    try std.testing.expectEqual(19, y3[0]);
    try std.testing.expectEqual(58, y3[2]);
    try std.testing.expectEqual(121, y3[4]);
    try std.testing.expectEqual(208, y3[6]);
    try std.testing.expectEqual(187, y3[8]);

    const x4 = try a.alloc(f64, 2 * n);
    defer a.free(x4);
    const y4 = try a.alloc(f64, 2 * n);
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

    sbmv(f64, .ColumnMajor, .Upper, n, k, alpha, A.ptr, k + 1, x4.ptr, -2, beta, y4.ptr, -2);

    try std.testing.expectEqual(19, y4[8]);
    try std.testing.expectEqual(58, y4[6]);
    try std.testing.expectEqual(121, y4[4]);
    try std.testing.expectEqual(208, y4[2]);
    try std.testing.expectEqual(187, y4[0]);

    const x5 = try a.alloc(f64, 2 * n);
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
        5,
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

    sbmv(f64, .RowMajor, .Lower, n, k, alpha, A.ptr, k + 1, x5.ptr, 2, beta, y5.ptr, 2);

    try std.testing.expectEqual(19, y5[0]);
    try std.testing.expectEqual(58, y5[2]);
    try std.testing.expectEqual(121, y5[4]);
    try std.testing.expectEqual(208, y5[6]);
    try std.testing.expectEqual(187, y5[8]);

    const x6 = try a.alloc(f64, 2 * n);
    defer a.free(x6);
    const y6 = try a.alloc(f64, 2 * n);
    defer a.free(y6);

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

    sbmv(f64, .RowMajor, .Lower, n, k, alpha, A.ptr, k + 1, x6.ptr, -2, beta, y6.ptr, -2);

    try std.testing.expectEqual(19, y6[8]);
    try std.testing.expectEqual(58, y6[6]);
    try std.testing.expectEqual(121, y6[4]);
    try std.testing.expectEqual(208, y6[2]);
    try std.testing.expectEqual(187, y6[0]);

    const x7 = try a.alloc(f64, 2 * n);
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
        5,
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

    sbmv(f64, .ColumnMajor, .Lower, n, k, alpha, A.ptr, k + 1, x7.ptr, 2, beta, y7.ptr, 2);

    try std.testing.expectEqual(13, y7[0]);
    try std.testing.expectEqual(46, y7[2]);
    try std.testing.expectEqual(103, y7[4]);
    try std.testing.expectEqual(184, y7[6]);
    try std.testing.expectEqual(169, y7[8]);

    const x8 = try a.alloc(f64, 2 * n);
    defer a.free(x8);
    const y8 = try a.alloc(f64, 2 * n);
    defer a.free(y8);

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

    sbmv(f64, .ColumnMajor, .Lower, n, k, alpha, A.ptr, k + 1, x8.ptr, -2, beta, y8.ptr, -2);

    try std.testing.expectEqual(13, y8[8]);
    try std.testing.expectEqual(46, y8[6]);
    try std.testing.expectEqual(103, y8[4]);
    try std.testing.expectEqual(184, y8[2]);
    try std.testing.expectEqual(169, y8[0]);
}
