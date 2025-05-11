const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const nrm2 = zml.linalg.blas.nrm2;

test nrm2 {
    const a = std.testing.allocator;

    const n = 1000;

    var x1 = try a.alloc(f64, n);
    defer a.free(x1);

    for (0..n) |i| {
        x1[i] = @floatFromInt(i + 1);
    }

    const result1 = nrm2(f64, n, x1.ptr, 1);
    try std.testing.expectApproxEqRel(18271.111077326415, result1, 0.0000000001);
    const result2 = nrm2(f64, n / 2, x1.ptr, 2);
    try std.testing.expectApproxEqRel(12909.9380323842, result2, 0.0000000001);

    var x2 = try a.alloc(cf64, n);
    defer a.free(x2);

    for (0..n) |i| {
        x2[i] = cf64.init(@floatFromInt(i + 1), @floatFromInt(-@as(isize, @intCast(i + 1))));
    }

    const result3 = nrm2(cf64, n, x2.ptr, 1);
    try std.testing.expectApproxEqRel(25839.253085180306, result3, 0.0000000001);
    const result4 = nrm2(cf64, n / 2, x2.ptr, 2);
    try std.testing.expectApproxEqRel(18257.409454793964, result4, 0.0000000001);
}
