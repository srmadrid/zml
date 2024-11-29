const std = @import("std");

// Numerical computing.
pub const ndarray = @import("ndarray/ndarray.zig");
pub const NDArray = ndarray.NDArray;
pub const BLAS = @import("ndarray/BLAS/BLAS.zig");
//pub const LAPACK = @import("ndarray/LAPACK/LAPACK.zig");

// Core symbolic constructs.
//pub const Expression = @import("expression/expression.zig").Expression;
//pub const Symbol = @import("symbol.zig").Symbol;
//pub const Element = @import("element.zig").Element;
//pub const Variable = @import("variable.zig").Variable;
//pub const Set = @import("set.zig").Set;

test {
    std.testing.refAllDeclsRecursive(@This());
}
