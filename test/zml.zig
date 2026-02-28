const std = @import("std");

const zml = @import("zml");

pub const float = @import("float.zig");
pub const cfloat = @import("cfloat.zig");

pub inline fn expectApproxEqAbs(expected: anytype, actual: anytype, tolerance: anytype) !void {
    if (expected != actual)
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
}

test {
    // Override test flags
    const test_all = false;

    // Individual test flags
    const test_types = true;
    const test_int = false;
    const test_float = false;
    const test_dyadic = false;
    const test_cfloat = false;
    const test_integer = false;
    const test_rational = false;
    const test_real = false;
    const test_complex = false;
    const test_constants = false;
    const test_numeric = false;
    const test_vector = false;
    const test_matrix = false;
    const test_array = false;
    const test_ops = false;
    const test_linalg = false;
    const test_autodiff = false;

    _ = test_int;
    _ = test_dyadic;
    _ = test_integer;
    _ = test_rational;
    _ = test_real;
    _ = test_complex;
    _ = test_constants;
    _ = test_numeric;
    _ = test_vector;
    _ = test_matrix;
    _ = test_array;
    _ = test_ops;
    _ = test_autodiff;

    if (test_all or test_types)
        _ = @import("types.zig");

    if (test_all or test_float)
        _ = @import("float.zig");

    if (test_all or test_cfloat)
        _ = @import("cfloat.zig");

    if (test_all or test_linalg)
        _ = @import("linalg.zig");
}
