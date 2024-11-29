const std = @import("std");

//
pub const getrf = @import("getrf.zig").getrf;

test {
    std.testing.refAllDeclsRecursive(@This());
}
