const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const iamax = zml.linalg.blas.iamax;

test iamax {
    const a: std.mem.Allocator = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = 0;
    }

    x1[127] = 1;
    x1[456] = 1;

    const result1 = iamax(n, x1.ptr, 1, .{}) catch unreachable;
    try std.testing.expectEqual(127, result1);
    // const result2 = iamax(n, x1.ptr, -1, .{}) catch unreachable;
    // try std.testing.expectEqual(0, result2);
    const result3 = iamax(n / 2, x1.ptr, 2, .{}) catch unreachable;
    try std.testing.expectEqual(228, result3);

    var x2 = try a.alloc(cf64, n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = cf64.init(0, 0);
    }

    x2[127] = cf64.init(1, 1);
    x2[456] = cf64.init(1, 1);

    const result4 = iamax(n, x2.ptr, 1, .{}) catch unreachable;
    try std.testing.expectEqual(127, result4);
    // const result5 = iamax(n, x2.ptr, -1, .{}) catch unreachable;
    // try std.testing.expectEqual(0, result5);
    const result6 = iamax(n / 2, x2.ptr, 2, .{}) catch unreachable;
    try std.testing.expectEqual(228, result6);
}
