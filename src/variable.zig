const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const Set = zml.Set;
const Expression = zml.Expression;

/// A mathematical variable.
pub const Variable = struct {
    /// Symbol of the variable.
    symbol: *Symbol,
    /// Set in which the variable lives.
    domain: *Set,
    /// Conditions that the elements satisfy.
    assumptions: Expression,
};
