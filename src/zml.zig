const std = @import("std");

pub const ndarray = @import("ndarray/ndarray.zig");
pub const NDArray = ndarray.NDArray;

test {
    std.testing.refAllDeclsRecursive(@This());
}
