const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const asum = zml.linalg.blas.asum;

test asum {
    const a = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
    }

    const result1 = asum(n, x1.ptr, 1, .{}) catch unreachable;
    try std.testing.expectEqual(500500, result1);
    const result2 = asum(n / 2, x1.ptr, 2, .{}) catch unreachable;
    try std.testing.expectEqual(250000, result2);

    var x2 = try a.alloc(cf64, n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
    }

    const result3 = asum(n, x2.ptr, 1, .{}) catch unreachable;
    try std.testing.expectEqual(1001000, result3);
    const result4 = asum(n / 2, x2.ptr, 2, .{}) catch unreachable;
    try std.testing.expectEqual(500000, result4);
}
