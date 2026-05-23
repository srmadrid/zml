const std = @import("std");

const zsl = @import("zsl");

pub const float = @import("float.zig");
pub const dyadic = @import("dyadic.zig");
pub const complex = @import("complex.zig");

pub const vector = @import("vector.zig");
pub const matrix = @import("matrix.zig");

pub inline fn expectApproxEqAbs(expected: anytype, actual: anytype, tolerance: anytype) !void {
    if (expected != actual)
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
}

pub fn randomNumber(comptime N: type, rand: std.Random) N {
    return zsl.numeric.cast(
        N,
        if (comptime zsl.meta.isComplex(N))
            zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
        else
            rand.float(f64),
    );
}

pub fn deinit(allocator: std.mem.Allocator, value: anytype) void {
    const V: type = zsl.meta.Child(@TypeOf(value));

    switch (comptime zsl.meta.domain(V)) {
        .numeric => {},
        else => value.deinit(allocator),
    }
}

test {
    // Override test flags
    const test_all = false;

    // Individual test flags
    const test_int = false;
    const test_float = false;
    const test_dyadic = false;
    const test_complex = false;
    const test_numeric = false;
    const test_vector = false;
    const test_matrix = false;
    const test_array = false;
    const test_ops = false;
    const test_linalg = true;
    const test_autodiff = false;

    _ = test_int;
    _ = test_numeric;
    _ = test_array;
    _ = test_ops;
    _ = test_autodiff;

    if (test_all or test_float)
        _ = @import("float.zig");

    if (test_all or test_dyadic)
        _ = @import("dyadic.zig");

    if (test_all or test_complex)
        _ = @import("complex.zig");

    if (test_all or test_vector)
        _ = @import("vector.zig");

    if (test_all or test_matrix)
        _ = @import("matrix.zig");

    if (test_all or test_linalg)
        _ = @import("linalg.zig");
}
