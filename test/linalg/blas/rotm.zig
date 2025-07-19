const std = @import("std");
const zml = @import("zml");
const rotm = zml.linalg.blas.rotm;

test rotm {
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

    const p1: []const f64 = &.{ 1, 2, 3, 4, 5 };
    rotm(n, x1.ptr, 1, x2.ptr, 1, p1.ptr, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + (i + 1))), x1[i]);
        try std.testing.expectEqual(@as(f64, @floatFromInt(5 * (i + 1) - (i + 1))), x2[i]);
        x1[i] = @floatFromInt(i + 1);
    }
    const p2: []const f64 = &.{ 0, 2, 3, 4, 5 };
    rotm(n, x1.ptr, 1, x3.ptr, -1, p2.ptr, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt((i + 1) + 4 * (i + 1))), x1[i]);
        try std.testing.expectEqual(@as(f64, @floatFromInt(4 * (n - i))), x3[i]);
        x1[i] = @floatFromInt(i + 1);
    }
    const p3: []const f64 = &.{ -1, 2, 3, 4, 5 };
    rotm(n / 2, x1.ptr, 2, x4.ptr, 2, p3.ptr, .{}) catch unreachable;
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 4 * (i + 1))), x1[i]);
            try std.testing.expectEqual(@as(f64, @floatFromInt(3 * (i + 1) + 5 * (i + 1))), x4[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x1[i]);
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x4[i]);
        }
    }
}
