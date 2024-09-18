const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const Expression = zml.Expression;

/// A
pub const Set = struct {
    elements: Expression,
    properties: []const Expression,
};
