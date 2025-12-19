const std = @import("std");

const zml = @import("zml");

pub const float = @import("float.zig");
pub const cfloat = @import("cfloat.zig");

pub inline fn expectApproxEqAbs(expected: anytype, actual: anytype, tolerance: anytype) !void {
    if (expected != actual)
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
}

test {
    const test_types = false;
    const test_float = true;
    const test_cfloat = true;
    const test_linalg = false;

    if (test_types) {
        _ = @import("types.zig");
    }

    if (test_float) {
        _ = @import("float.zig");
    }

    if (test_cfloat) {
        _ = @import("cfloat.zig");
    }

    if (test_linalg) {
        _ = @import("linalg.zig");
    }
}
