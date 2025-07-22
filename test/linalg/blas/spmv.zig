const std = @import("std");
const zml = @import("zml");
const spmv = zml.linalg.blas.spmv;

test spmv {
    const a = std.testing.allocator;

    const n = 5;
    const alpha = 2;
    const beta = 3;

    const A = try a.alloc(f64, n * n);
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

    spmv(.row_major, .upper, n, alpha, A.ptr, x1.ptr, 2, beta, y1.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(113, y1[0]);
    try std.testing.expectEqual(230, y1[2]);
    try std.testing.expectEqual(311, y1[4]);
    try std.testing.expectEqual(362, y1[6]);
    try std.testing.expectEqual(395, y1[8]);

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

    spmv(.row_major, .upper, n, alpha, A.ptr, x2.ptr, -2, beta, y2.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(113, y2[8]);
    try std.testing.expectEqual(230, y2[6]);
    try std.testing.expectEqual(311, y2[4]);
    try std.testing.expectEqual(362, y2[2]);
    try std.testing.expectEqual(395, y2[0]);

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

    spmv(.col_major, .upper, n, alpha, A.ptr, x3.ptr, 2, beta, y3.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(203, y3[0]);
    try std.testing.expectEqual(236, y3[2]);
    try std.testing.expectEqual(275, y3[4]);
    try std.testing.expectEqual(332, y3[6]);
    try std.testing.expectEqual(425, y3[8]);

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

    spmv(.col_major, .upper, n, alpha, A.ptr, x4.ptr, -2, beta, y4.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(203, y4[8]);
    try std.testing.expectEqual(236, y4[6]);
    try std.testing.expectEqual(275, y4[4]);
    try std.testing.expectEqual(332, y4[2]);
    try std.testing.expectEqual(425, y4[0]);

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

    spmv(.row_major, .lower, n, alpha, A.ptr, x5.ptr, 2, beta, y5.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(203, y5[0]);
    try std.testing.expectEqual(236, y5[2]);
    try std.testing.expectEqual(275, y5[4]);
    try std.testing.expectEqual(332, y5[6]);
    try std.testing.expectEqual(425, y5[8]);

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

    spmv(.row_major, .lower, n, alpha, A.ptr, x6.ptr, -2, beta, y6.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(203, y6[8]);
    try std.testing.expectEqual(236, y6[6]);
    try std.testing.expectEqual(275, y6[4]);
    try std.testing.expectEqual(332, y6[2]);
    try std.testing.expectEqual(425, y6[0]);

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

    spmv(.col_major, .lower, n, alpha, A.ptr, x7.ptr, 2, beta, y7.ptr, 2, .{}) catch unreachable;

    try std.testing.expectEqual(113, y7[0]);
    try std.testing.expectEqual(230, y7[2]);
    try std.testing.expectEqual(311, y7[4]);
    try std.testing.expectEqual(362, y7[6]);
    try std.testing.expectEqual(395, y7[8]);

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

    spmv(.col_major, .lower, n, alpha, A.ptr, x8.ptr, -2, beta, y8.ptr, -2, .{}) catch unreachable;

    try std.testing.expectEqual(113, y8[8]);
    try std.testing.expectEqual(230, y8[6]);
    try std.testing.expectEqual(311, y8[4]);
    try std.testing.expectEqual(362, y8[2]);
    try std.testing.expectEqual(395, y8[0]);
}
