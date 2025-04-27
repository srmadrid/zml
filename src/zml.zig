const std = @import("std");

// Core types and operations.
pub const core = @import("core.zig");

// Numerical computing.
//pub const ndarray = @import("ndarray/ndarray.zig");
//pub const NDArray = ndarray.NDArray;
//pub const linalg = @import("linalg.zig");

// Core symbolic constructs.
//pub const Expression = @import("expression/expression.zig").Expression;
//pub const Symbol = @import("symbol.zig").Symbol;
//pub const Element = @import("element.zig").Element;
//pub const Variable = @import("variable.zig").Variable;
//pub const Set = @import("set.zig").Set;

test {
    std.testing.refAllDeclsRecursive(@This());
}
