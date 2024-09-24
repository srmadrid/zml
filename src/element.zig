const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const Set = zml.Set;
const Expression = zml.Expression;

/// A mathematical element.
pub const Element = struct {
    /// Symbol of the variable.
    symbol: *Symbol,
    /// Set in which the element lives.
    domain: Set,
    /// The element itself. NOT FINAL. MAYBE HAVE THE NAME OF THE SYMBOL BE THE ELEMENT
    element: *anyopaque,
};
