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
pub const exp10 = ops.exp10;
pub const exp10_ = ops.exp10_;
pub const exp10_to = ops.exp10_to;
pub const exp2 = ops.exp2;
pub const exp2_ = ops.exp2_;
pub const exp2_to = ops.exp2_to;
pub const exp10m1 = ops.exp10m1;
pub const exp10m1_ = ops.exp10m1_;
pub const exp10m1_to = ops.exp10m1_to;
pub const exp2m1 = ops.exp2m1;
pub const exp2m1_ = ops.exp2m1_;
pub const exp2m1_to = ops.exp2m1_to;
pub const expm1 = ops.expm1;
pub const expm1_ = ops.expm1_;
pub const expm1_to = ops.expm1_to;
pub const log = ops.log;
pub const log_ = ops.log_;
pub const log_to = ops.log_to;
pub const log10 = ops.log10;
pub const log10_ = ops.log10_;
pub const log10_to = ops.log10_to;
pub const log2 = ops.log2;
pub const log2_ = ops.log2_;
pub const log2_to = ops.log2_to;
pub const log10p1 = ops.log10p1;
pub const log10p1_ = ops.log10p1_;
pub const log10p1_to = ops.log10p1_to;
pub const log2p1 = ops.log2p1;
pub const log2p1_ = ops.log2p1_;
pub const log2p1_to = ops.log2p1_to;
pub const log1p = ops.log1p;
pub const log1p_ = ops.log1p_;
pub const log1p_to = ops.log1p_to;

// Power functions
pub const pow = ops.pow;
pub const pow_ = ops.pow_;
pub const pow_to = ops.pow_to;
pub const sqrt = ops.sqrt;
pub const sqrt_ = ops.sqrt_;
pub const sqrt_to = ops.sqrt_to;
pub const cbrt = ops.cbrt;
pub const cbrt_ = ops.cbrt_;
pub const cbrt_to = ops.cbrt_to;
pub const hypot = ops.hypot;
pub const hypot_ = ops.hypot_;
pub const hypot_to = ops.hypot_to;

// Trigonometric functions
pub const sin = ops.sin;
pub const sin_ = ops.sin_;
pub const sin_to = ops.sin_to;
pub const cos = ops.cos;
pub const cos_ = ops.cos_;
pub const cos_to = ops.cos_to;
pub const tan = ops.tan;
pub const tan_ = ops.tan_;
pub const tan_to = ops.tan_to;
pub const asin = ops.asin;
pub const asin_ = ops.asin_;
pub const asin_to = ops.asin_to;
pub const acos = ops.acos;
pub const acos_ = ops.acos_;
pub const acos_to = ops.acos_to;
pub const atan = ops.atan;
pub const atan_ = ops.atan_;
pub const atan_to = ops.atan_to;
pub const atan2 = ops.atan2;
pub const atan2_ = ops.atan2_;
pub const atan2_to = ops.atan2_to;
pub const sinpi = ops.sinpi;
pub const sinpi_ = ops.sinpi_;
pub const sinpi_to = ops.sinpi_to;
pub const cospi = ops.cospi;
pub const cospi_ = ops.cospi_;
pub const cospi_to = ops.cospi_to;
pub const tanpi = ops.tanpi;
pub const tanpi_ = ops.tanpi_;
pub const tanpi_to = ops.tanpi_to;
pub const asinpi = ops.asinpi;
pub const asinpi_ = ops.asinpi_;
pub const asinpi_to = ops.asinpi_to;
pub const acospi = ops.acospi;
pub const acospi_ = ops.acospi_;
pub const acospi_to = ops.acospi_to;
pub const atanpi = ops.atanpi;
pub const atanpi_ = ops.atanpi_;
pub const atanpi_to = ops.atanpi_to;
pub const atan2pi = ops.atan2pi;
pub const atan2pi_ = ops.atan2pi_;
pub const atan2pi_to = ops.atan2pi_to;

// Hyperbolic functions
pub const sinh = ops.sinh;
pub const sinh_ = ops.sinh_;
pub const sinh_to = ops.sinh_to;
pub const cosh = ops.cosh;
pub const cosh_ = ops.cosh_;
pub const cosh_to = ops.cosh_to;
pub const tanh = ops.tanh;
pub const tanh_ = ops.tanh_;
pub const tanh_to = ops.tanh_to;
pub const asinh = ops.asinh;
pub const asinh_ = ops.asinh_;
pub const asinh_to = ops.asinh_to;
pub const acosh = ops.acosh;
pub const acosh_ = ops.acosh_;
pub const acosh_to = ops.acosh_to;
pub const atanh = ops.atanh;
pub const atanh_ = ops.atanh_;
pub const atanh_to = ops.atanh_to;

// Error and gamma functions
pub const erf = ops.erf;
pub const erf_ = ops.erf_;
pub const erf_to = ops.erf_to;
pub const erfc = ops.erfc;
pub const erfc_ = ops.erfc_;
pub const erfc_to = ops.erfc_to;
pub const gamma = ops.gamma;
pub const gamma_ = ops.gamma_;
pub const gamma_to = ops.gamma_to;
pub const lgamma = ops.lgamma;
pub const lgamma_ = ops.lgamma_;
pub const lgamma_to = ops.lgamma_to;

// Nearest integer operations
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
