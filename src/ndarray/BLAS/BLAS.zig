const std = @import("std");

pub const asum = @import("asum.zig").asum;

test {
    std.testing.refAllDeclsRecursive(@This());
}
