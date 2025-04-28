const std = @import("std");

pub const types = @import("core/types.zig");
pub const math = @import("core/math.zig");
pub const ops = @import("core/ops.zig");

test {
    const test_types = true;
    const test_math = true;

    if (test_types) {
        _ = @import("core/types.zig");
    }

    if (test_math) {
        _ = @import("core/math.zig");
    }
}
