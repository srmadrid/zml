const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const copy = zml.linalg.blas.copy;

test copy {
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
        x2[i] = 0;
        x3[i] = 0;
        x4[i] = 0;
    }

    copy(n, x1.ptr, 1, x2.ptr, 1, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x2[i]);
    }
    copy(n, x1.ptr, 1, x3.ptr, -1, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(n - i)), x3[i]);
    }
    copy(n / 2, x1.ptr, 2, x4.ptr, 2, .{}) catch unreachable;
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x4[i]);
        } else {
            try std.testing.expectEqual(0, x4[i]);
        }
    }

    var x5 = try a.alloc(cf64, n);
    defer a.free(x5);
    var x6 = try a.alloc(cf64, n);
    defer a.free(x6);
    var x7 = try a.alloc(cf64, n);
    defer a.free(x7);
    var x8 = try a.alloc(cf64, n);
    defer a.free(x8);

    for (0..n) |i| {
        x5[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
        x6[i] = cf64.init(0, 0);
        x7[i] = cf64.init(0, 0);
        x8[i] = cf64.init(0, 0);
    }

    copy(n, x5.ptr, 1, x6.ptr, 1, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x6[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x6[i].im);
    }
    copy(n, x5.ptr, 1, x7.ptr, -1, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(n - i)), x7[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(n - i)))), x7[i].im);
    }
    copy(n / 2, x5.ptr, 2, x8.ptr, 2, .{}) catch unreachable;
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x8[i].im);
        } else {
            try std.testing.expectEqual(0, x8[i].re);
            try std.testing.expectEqual(0, x8[i].im);
        }
    }
}
