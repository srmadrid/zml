const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const axpy = zml.linalg.blas.axpy;

test axpy {
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

    axpy(f64, n, 2, x1.ptr, 1, x2.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + (i + 1))), x2[i]);
    }
    axpy(f64, n, 2, x1.ptr, 1, x3.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i) + (n - i))), x3[i]);
    }
    axpy(f64, n / 2, 2, x1.ptr, 2, x4.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + (i + 1))), x4[i]);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x4[i]);
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
        x5[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
        x6[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
        x7[i] = cf64.init(@floatFromInt(n - i), @floatFromInt(-@as(isize, @intCast((n - i)))));
        x8[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast((i + 1)))));
    }

    axpy(cf64, n, cf64.init(2, 3), x5.ptr, 1, x6.ptr, 1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1) + (i + 1))), x6[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)) - @as(isize, @intCast(i + 1)))), x6[i].im);
    }
    axpy(cf64, n, cf64.init(2, 3), x5.ptr, 1, x7.ptr, -1);
    for (0..n) |i| {
        try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (n - i) + 3 * (n - i) + (n - i))), x7[i].re);
        try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(n - i)) + 3 * @as(isize, @intCast(n - i)) - @as(isize, @intCast(n - i)))), x7[i].im);
    }
    axpy(cf64, n / 2, cf64.init(2, 3), x5.ptr, 2, x8.ptr, 2);
    for (0..n) |i| {
        if (i % 2 == 0) {
            try std.testing.expectEqual(@as(f64, @floatFromInt(2 * (i + 1) + 3 * (i + 1) + (i + 1))), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-2 * @as(isize, @intCast(i + 1)) + 3 * @as(isize, @intCast(i + 1)) - @as(isize, @intCast(i + 1)))), x8[i].im);
        } else {
            try std.testing.expectEqual(@as(f64, @floatFromInt(i + 1)), x8[i].re);
            try std.testing.expectEqual(@as(f64, @floatFromInt(-@as(isize, @intCast(i + 1)))), x8[i].im);
        }
    }
}
