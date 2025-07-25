const std = @import("std");
const zml = @import("zml");
const ger = zml.linalg.blas.ger;

test ger {
    const a = std.testing.allocator;

    const m = 4;
    const n = 5;
    const alpha: f64 = 2;

    const A = try a.alloc(f64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(f64, 2 * m);
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

    ger(.row_major, m, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(3, A[0]);
    try std.testing.expectEqual(6, A[1]);
    try std.testing.expectEqual(9, A[2]);
    try std.testing.expectEqual(12, A[3]);
    try std.testing.expectEqual(15, A[4]);
    try std.testing.expectEqual(10, A[5]);
    try std.testing.expectEqual(15, A[6]);
    try std.testing.expectEqual(20, A[7]);
    try std.testing.expectEqual(25, A[8]);
    try std.testing.expectEqual(30, A[9]);
    try std.testing.expectEqual(17, A[10]);
    try std.testing.expectEqual(24, A[11]);
    try std.testing.expectEqual(31, A[12]);
    try std.testing.expectEqual(38, A[13]);
    try std.testing.expectEqual(45, A[14]);
    try std.testing.expectEqual(24, A[15]);
    try std.testing.expectEqual(33, A[16]);
    try std.testing.expectEqual(42, A[17]);
    try std.testing.expectEqual(51, A[18]);
    try std.testing.expectEqual(60, A[19]);

    const x2 = try a.alloc(f64, 2 * m);
    defer a.free(x2);
    const y2 = try a.alloc(f64, 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]f64{
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

    ger(.row_major, m, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(5, A[0]);
    try std.testing.expectEqual(10, A[1]);
    try std.testing.expectEqual(15, A[2]);
    try std.testing.expectEqual(20, A[3]);
    try std.testing.expectEqual(25, A[4]);
    try std.testing.expectEqual(14, A[5]);
    try std.testing.expectEqual(23, A[6]);
    try std.testing.expectEqual(32, A[7]);
    try std.testing.expectEqual(41, A[8]);
    try std.testing.expectEqual(50, A[9]);
    try std.testing.expectEqual(23, A[10]);
    try std.testing.expectEqual(36, A[11]);
    try std.testing.expectEqual(49, A[12]);
    try std.testing.expectEqual(62, A[13]);
    try std.testing.expectEqual(75, A[14]);
    try std.testing.expectEqual(32, A[15]);
    try std.testing.expectEqual(49, A[16]);
    try std.testing.expectEqual(66, A[17]);
    try std.testing.expectEqual(83, A[18]);
    try std.testing.expectEqual(100, A[19]);

    const x3 = try a.alloc(f64, 2 * m);
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

    ger(.col_major, m, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(7, A[0]);
    try std.testing.expectEqual(14, A[1]);
    try std.testing.expectEqual(21, A[2]);
    try std.testing.expectEqual(28, A[3]);
    try std.testing.expectEqual(29, A[4]);
    try std.testing.expectEqual(22, A[5]);
    try std.testing.expectEqual(35, A[6]);
    try std.testing.expectEqual(48, A[7]);
    try std.testing.expectEqual(47, A[8]);
    try std.testing.expectEqual(62, A[9]);
    try std.testing.expectEqual(41, A[10]);
    try std.testing.expectEqual(60, A[11]);
    try std.testing.expectEqual(57, A[12]);
    try std.testing.expectEqual(78, A[13]);
    try std.testing.expectEqual(99, A[14]);
    try std.testing.expectEqual(64, A[15]);
    try std.testing.expectEqual(59, A[16]);
    try std.testing.expectEqual(86, A[17]);
    try std.testing.expectEqual(113, A[18]);
    try std.testing.expectEqual(140, A[19]);

    const x4 = try a.alloc(f64, 2 * m);
    defer a.free(x4);
    const y4 = try a.alloc(f64, 2 * n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]f64{
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

    ger(.col_major, m, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(9, A[0]);
    try std.testing.expectEqual(18, A[1]);
    try std.testing.expectEqual(27, A[2]);
    try std.testing.expectEqual(36, A[3]);
    try std.testing.expectEqual(33, A[4]);
    try std.testing.expectEqual(30, A[5]);
    try std.testing.expectEqual(47, A[6]);
    try std.testing.expectEqual(64, A[7]);
    try std.testing.expectEqual(53, A[8]);
    try std.testing.expectEqual(74, A[9]);
    try std.testing.expectEqual(59, A[10]);
    try std.testing.expectEqual(84, A[11]);
    try std.testing.expectEqual(65, A[12]);
    try std.testing.expectEqual(94, A[13]);
    try std.testing.expectEqual(123, A[14]);
    try std.testing.expectEqual(96, A[15]);
    try std.testing.expectEqual(69, A[16]);
    try std.testing.expectEqual(106, A[17]);
    try std.testing.expectEqual(143, A[18]);
    try std.testing.expectEqual(180, A[19]);
}
