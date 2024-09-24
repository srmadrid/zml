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
    /// THIS CAN BE A POINTER OR NOT: IT BEING A POINTER IS EASY; IT NOT BEING
    /// A POINTER IS MORE COMPLEX BUT IS FINE, SINCE IT IS STILL STORED AS A
    /// DEPENDENT IN ITS SYMBOL AND THERE IS ALSO A POINTER TO THE SET'S SYMBOL
    domain: Set,
    /// Conditions that the elements satisfy.
    assumptions: Expression,
};
