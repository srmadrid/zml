const std = @import("std");

pub const types = @import("types.zig");
pub const Order = types.Order;
pub const scast = types.scast;
pub const cast = types.cast;

pub const int = @import("int.zig");
pub const float = @import("float.zig");
pub const cfloat = @import("cfloat.zig");
pub const cf16 = cfloat.cf16;
pub const cf32 = cfloat.cf32;
pub const cf64 = cfloat.cf64;
pub const cf80 = cfloat.cf80;
pub const cf128 = cfloat.cf128;
pub const comptime_complex = cfloat.comptime_complex;
pub const integer = @import("integer.zig");
pub const Integer = integer.Integer;
pub const rational = @import("rational.zig");
pub const Rational = rational.Rational;
pub const real = @import("real.zig");
pub const Real = real.Real;
pub const complex = @import("complex.zig");
pub const Complex = complex.Complex;
pub const array = @import("array.zig");
pub const Array = array.Array;

const constants = @import("constants.zig");
pub const pi = constants.pi;

const ops = @import("ops.zig");
pub const add = ops.add;
pub const add_ = ops.add_;
pub const add_to = ops.add_to;
pub const sub = ops.sub;
pub const sub_ = ops.sub_;
pub const sub_to = ops.sub_to;
pub const mul = ops.mul;
pub const mul_ = ops.mul_;
pub const mul_to = ops.mul_to;
pub const div = ops.div;
pub const div_ = ops.div_;
pub const div_to = ops.div_to;

pub const eq = ops.eq;
pub const ne = ops.ne;
pub const lt = ops.lt;
pub const le = ops.le;
pub const gt = ops.gt;
pub const ge = ops.ge;

pub const max = ops.max;
pub const min = ops.min;

// Basic operations
pub const abs = ops.abs;
pub const abs_ = ops.abs_;
pub const abs_to = ops.abs_to;

// Exponential functions
pub const exp = ops.exp;
pub const exp_ = ops.exp_;
pub const exp_to = ops.exp_to;

pub const ceil = ops.ceil;
pub const ceil_ = ops.ceil_;
pub const ceil_to = ops.ceil_to;

pub const linalg = @import("linalg.zig");

// Symbolic system.
//pub const Expression = @import("expression/expression.zig").Expression;
//pub const Symbol = @import("symbol.zig").Symbol;
//pub const Element = @import("element.zig").Element;
//pub const Variable = @import("variable.zig").Variable;
//pub const Set = @import("set.zig").Set;
