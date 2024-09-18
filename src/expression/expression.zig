const std = @import("std");
const zml = @import("../zml.zig");
const ExpressionTree = @import("exprtree.zig").ExpressionTree;

/// A
pub const Expression = struct {
    string: []const u8,
    tree: ExpressionTree,
};

test {
    _ = @import("exprtree.zig");
    std.testing.refAllDeclsRecursive(@This());
}
