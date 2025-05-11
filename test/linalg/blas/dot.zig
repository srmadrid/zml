const std = @import("std");
const zml = @import("zml");
const dot = zml.linalg.blas.dot;

test dot {
    const a: std.mem.Allocator = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);
    var x2 = try a.alloc(f64, n);
    defer a.free(x2);
    var x3 = try a.alloc(f64, n);
    defer a.free(x3);
    var x4 = try a.alloc(f64, n);
    defer a.free(x4);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
        x2[i] = @floatFromInt(i + 1);
        x3[i] = @floatFromInt(n - i);
        x4[i] = @floatFromInt(i + 1);
    }

    const result1 = dot(f64, n, x1.ptr, 1, x2.ptr, 1);
    try std.testing.expectEqual(333833500, result1);
    const result2 = dot(f64, n, x1.ptr, 1, x3.ptr, -1);
    try std.testing.expectEqual(333833500, result2);
    const result3 = dot(f64, n / 2, x1.ptr, 2, x4.ptr, 2);
    try std.testing.expectEqual(166666500, result3);
}
