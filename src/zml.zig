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
pub const sub = ops.sub;
pub const sub_ = ops.sub_;
pub const mul = ops.mul;
pub const mul_ = ops.mul_;
pub const div = ops.div;
pub const div_ = ops.div_;

pub const eq = ops.eq;
pub const eq_ = ops.eq_;
pub const ne = ops.ne;
pub const ne_ = ops.ne_;
pub const lt = ops.lt;
pub const lt_ = ops.lt_;
pub const le = ops.le;
pub const le_ = ops.le_;
pub const gt = ops.gt;
pub const gt_ = ops.gt_;
pub const ge = ops.ge;
pub const ge_ = ops.ge_;

pub const max = ops.max;
pub const min = ops.min;

// Basic operations
pub const abs = ops.abs;
pub const abs_ = ops.abs_;

// Exponential functions
pub const exp = ops.exp;
pub const exp_ = ops.exp_;
pub const exp10 = ops.exp10;
pub const exp10_ = ops.exp10_;
pub const exp2 = ops.exp2;
pub const exp2_ = ops.exp2_;
pub const exp10m1 = ops.exp10m1;
pub const exp10m1_ = ops.exp10m1_;
pub const exp2m1 = ops.exp2m1;
pub const exp2m1_ = ops.exp2m1_;
pub const expm1 = ops.expm1;
pub const expm1_ = ops.expm1_;
pub const log = ops.log;
pub const log_ = ops.log_;
pub const log10 = ops.log10;
pub const log10_ = ops.log10_;
pub const log2 = ops.log2;
pub const log2_ = ops.log2_;
pub const log10p1 = ops.log10p1;
pub const log10p1_ = ops.log10p1_;
pub const log2p1 = ops.log2p1;
pub const log2p1_ = ops.log2p1_;
pub const log1p = ops.log1p;
pub const log1p_ = ops.log1p_;

// Power functions
pub const pow = ops.pow;
pub const pow_ = ops.pow_;
pub const sqrt = ops.sqrt;
pub const sqrt_ = ops.sqrt_;
pub const cbrt = ops.cbrt;
pub const cbrt_ = ops.cbrt_;
pub const hypot = ops.hypot;
pub const hypot_ = ops.hypot_;

// Trigonometric functions
pub const sin = ops.sin;
pub const sin_ = ops.sin_;
pub const cos = ops.cos;
pub const cos_ = ops.cos_;
pub const tan = ops.tan;
pub const tan_ = ops.tan_;
pub const asin = ops.asin;
pub const asin_ = ops.asin_;
pub const acos = ops.acos;
pub const acos_ = ops.acos_;
pub const atan = ops.atan;
pub const atan_ = ops.atan_;
pub const atan2 = ops.atan2;
pub const atan2_ = ops.atan2_;
pub const sinpi = ops.sinpi;
pub const sinpi_ = ops.sinpi_;
pub const cospi = ops.cospi;
pub const cospi_ = ops.cospi_;
pub const tanpi = ops.tanpi;
pub const tanpi_ = ops.tanpi_;
pub const asinpi = ops.asinpi;
pub const asinpi_ = ops.asinpi_;
pub const acospi = ops.acospi;
pub const acospi_ = ops.acospi_;
pub const atanpi = ops.atanpi;
pub const atanpi_ = ops.atanpi_;
pub const atan2pi = ops.atan2pi;
pub const atan2pi_ = ops.atan2pi_;

// Hyperbolic functions
pub const sinh = ops.sinh;
pub const sinh_ = ops.sinh_;
pub const cosh = ops.cosh;
pub const cosh_ = ops.cosh_;
pub const tanh = ops.tanh;
pub const tanh_ = ops.tanh_;
pub const asinh = ops.asinh;
pub const asinh_ = ops.asinh_;
pub const acosh = ops.acosh;
pub const acosh_ = ops.acosh_;
pub const atanh = ops.atanh;
pub const atanh_ = ops.atanh_;

// Error and gamma functions
pub const erf = ops.erf;
pub const erf_ = ops.erf_;
pub const erfc = ops.erfc;
pub const erfc_ = ops.erfc_;
pub const gamma = ops.gamma;
pub const gamma_ = ops.gamma_;
pub const lgamma = ops.lgamma;
pub const lgamma_ = ops.lgamma_;

// Nearest integer operations
pub const ceil = ops.ceil;
pub const ceil_ = ops.ceil_;

pub const set = ops.set;

pub const linalg = @import("linalg.zig");

// Symbolic system.
//pub const Expression = @import("expression/expression.zig").Expression;
//pub const Symbol = @import("symbol.zig").Symbol;
//pub const Element = @import("element.zig").Element;
//pub const Variable = @import("variable.zig").Variable;
//pub const Set = @import("set.zig").Set;
