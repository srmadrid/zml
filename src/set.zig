const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const Expression = zml.Expression;

/// A mathematical set.
pub const Set = struct {
    /// Symbol of the set.
    symbol: *Symbol,
    /// Variable from which the set is created. The domain of the set is
    /// contained here, and if it is equal to the set itself, and no condition
    /// is passed, then the set is only declared, i.e., it is equivalent to
    /// just saying ``Let $G$ be a set'' without explicitly defining $G$ itself.
    variable: Symbol,
    /// Conditions that the elements satisfy.
    condition: Expression,

    /// Initialization functions
    pub const init = struct {};
};
