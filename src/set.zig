const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const Expression = zml.Expression;

/// A mathematical set.
pub const Set = struct {
    /// Symbol of the set.
    symbol: Symbol,
    variable: Symbol,
    domain: Symbol,
    condition: Expression,

    /// Initialization functions
    pub const init = struct {};
};
