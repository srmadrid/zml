const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const scal = zml.linalg.blas.scal;

test scal {
    const a: std.mem.Allocator = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);
    var x2 = try a.alloc(f64, n);
    defer a.free(x2);
    var x3 = try a.alloc(f64, n);
    defer a.free(x3);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
        x2[i] = @floatFromInt(n - i);
        x3[i] = @floatFromInt(i + 1);
    }

    scal(n, 2, x1.ptr, 1, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1))), x1[i]);
    }
    // scal(n, 2, x2.ptr, -1, .{}) catch unreachable;
    // for (0..n) |i| {
    //     try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i))), x2[i]);
    // }
    scal(n / 2, 2, x3.ptr, 2, .{}) catch unreachable;
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1))), x3[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x3[i]);
        }
    }

    var x4 = try a.alloc(cf64, n);
    defer a.free(x4);
    var x5 = try a.alloc(cf64, n);
    defer a.free(x5);
    var x6 = try a.alloc(cf64, n);
    defer a.free(x6);

    for (0..n) |i| {
        x4[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
        x5[i] = cf64.init(@floatFromInt(n - i), @floatFromInt(-@as(isize, @intCast((n - i)))));
        x6[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
    }

    scal(n, cf64.init(2, 3), x4.ptr, 1, .{}) catch unreachable;
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x4[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)))), x4[i].im);
    }
    // scal(n, cf64.init(2, 3), x5.ptr, -1, .{}) catch unreachable;
    // for (0..n) |i| {
    //     try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i) + 3 * (n - i))), x5[i].re);
    //     try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(n - i)) + 3 * @as(isize, @intCast(n - i)))), x5[i].im);
    // }
    scal(n / 2, cf64.init(2, 3), x6.ptr, 2, .{}) catch unreachable;
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1))), x6[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)))), x6[i].im);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x6[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x6[i].im);
        }
    }
}
