const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const iamin = zml.linalg.blas.iamin;

test iamin {
    const a: std.mem.Allocator = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = std.math.floatMax(f64);
    }

    x1[127] = 0;
    x1[456] = 0;

    const result1 = iamin(f64, n, x1.ptr, 1);
    try std.testing.expectEqual(127, result1);
    const result2 = iamin(f64, n, x1.ptr, -1);
    try std.testing.expectEqual(0, result2);
    const result3 = iamin(f64, n / 2, x1.ptr, 2);
    try std.testing.expectEqual(228, result3);

    var x2 = try a.alloc(cf64, n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = cf64.init(std.math.floatMax(f64), std.math.floatMax(f64));
    }

    x2[127] = cf64.init(0, 0);
    x2[456] = cf64.init(0, 0);

    const result4 = iamin(cf64, n, x2.ptr, 1);
    try std.testing.expectEqual(127, result4);
    const result5 = iamin(cf64, n, x2.ptr, -1);
    try std.testing.expectEqual(0, result5);
    const result6 = iamin(cf64, n / 2, x2.ptr, 2);
    try std.testing.expectEqual(228, result6);
}
