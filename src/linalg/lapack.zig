const std = @import("std");

//
pub const getrf = @import("lapack/getrf.zig").getrf;

test {
    std.testing.refAllDeclsRecursive(@This());
}
