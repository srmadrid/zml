const std = @import("std");

// Core types and operations.
pub const types = @import("types.zig");
pub const ops = @import("ops.zig");
pub const float = @import("float.zig");
pub const cfloat = @import("cfloat.zig");
pub const cf16 = cfloat.cf16;
pub const cf32 = cfloat.cf32;
pub const cf64 = cfloat.cf64;
pub const cf80 = cfloat.cf80;
pub const cf128 = cfloat.cf128;
pub const integer = @import("integer.zig");
pub const Integer = integer.Integer;
pub const IntegerManaged = integer.IntegerManaged;
pub const rational = @import("rational.zig");
pub const Rational = rational.Rational;
pub const RationalManaged = rational.RationalManaged;
pub const real = @import("real.zig");
pub const Real = real.Real;
pub const RealManaged = real.RealManaged;
pub const complex = @import("complex.zig");
pub const Complex = complex.Complex;
pub const ComplexManaged = complex.ComplexManaged;

// Numerical computing.
//pub const ndarray = @import("ndarray/ndarray.zig");
//pub const NDArray = ndarray.NDArray;
pub const linalg = @import("linalg.zig");

// Core symbolic constructs.
//pub const Expression = @import("expression/expression.zig").Expression;
//pub const Symbol = @import("symbol.zig").Symbol;
//pub const Element = @import("element.zig").Element;
//pub const Variable = @import("variable.zig").Variable;
//pub const Set = @import("set.zig").Set;
