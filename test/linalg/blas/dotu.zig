const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const dotu = zml.linalg.blas.dotu;

test dotu {
    const a: std.mem.Allocator = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(cf64, n);
    defer a.free(x1);
    var x2 = try a.alloc(cf64, n);
    defer a.free(x2);
    var x3 = try a.alloc(cf64, n);
    defer a.free(x3);
    var x4 = try a.alloc(cf64, n);
    defer a.free(x4);

    for (0..n) |i| {
        x1[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
        x2[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
        x3[i] = cf64.init(@floatFromInt(n - i), @floatFromInt(-@as(isize, @intCast(n - i))));
        x4[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
    }

    const result1 = dotu(n, x1.ptr, 1, x2.ptr, 1, .{}) catch unreachable;
    try std.testing.expectEqual(0, result1.re);
    try std.testing.expectEqual(-667667000, result1.im);
    const result2 = dotu(n, x1.ptr, 1, x3.ptr, -1, .{}) catch unreachable;
    try std.testing.expectEqual(0, result2.re);
    try std.testing.expectEqual(-667667000, result2.im);
    const result3 = dotu(n / 2, x1.ptr, 2, x4.ptr, 2, .{}) catch unreachable;
    try std.testing.expectEqual(0, result3.re);
    try std.testing.expectEqual(-333333000, result3.im);
}
