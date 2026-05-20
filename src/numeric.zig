//! Namespace for numeric types and operations.

const constants = @import("numeric/constants.zig");
const ops = @import("numeric/ops.zig");

// Constants
pub const zero = constants.zero;
pub const one = constants.one;
pub const two = constants.two;

// Utilities
pub const Coerce = ops.Coerce;
pub const cast = ops.cast;
pub const set = ops.set;

// Basic operations
pub const Abs = ops.Abs;
pub const abs = ops.abs;
pub const abs_ = ops.abs_;
pub const Abs1 = ops.Abs1;
pub const abs1 = ops.abs1;
pub const abs1_ = ops.abs1_;
pub const Abs2 = ops.Abs2;
pub const abs2 = ops.abs2;
pub const abs2_ = ops.abs2_;
pub const Neg = ops.Neg;
pub const neg = ops.neg;
pub const neg_ = ops.neg_;
pub const Re = ops.Re;
pub const re = ops.re;
pub const Im = ops.Im;
pub const im = ops.im;
pub const Conj = ops.Conj;
pub const conj = ops.conj;
pub const conj_ = ops.conj_;
pub const Sign = ops.Sign;
pub const sign = ops.sign;
// pub const copysign = numops.copysign;

// Arithmetic operations
pub const Add = ops.Add;
pub const add = ops.add;
pub const add_ = ops.add_;
pub const atomicAdd_ = ops.atomicAdd_;
pub const Sub = ops.Sub;
pub const sub = ops.sub;
pub const sub_ = ops.sub_;
pub const Mul = ops.Mul;
pub const mul = ops.mul;
pub const mul_ = ops.mul_;
pub const Fma = ops.Fma;
pub const fma = ops.fma;
pub const fma_ = ops.fma_;
pub const Div = ops.Div;
pub const div = ops.div;
pub const div_ = ops.div_;

// Comparison operations
// pub const cmp = ops.cmp;
pub const eq = ops.eq;
pub const ne = ops.ne;
pub const lt = ops.lt;
pub const le = ops.le;
pub const gt = ops.gt;
pub const ge = ops.ge;
pub const Max = ops.Max;
pub const max = ops.max;
pub const max_ = ops.max_;
pub const Min = ops.Min;
pub const min = ops.min;
pub const min_ = ops.min_;

// Exponential functions
pub const Exp = ops.Exp;
pub const exp = ops.exp;
pub const exp_ = ops.exp_;
pub const Ln = ops.Ln;
pub const ln = ops.ln;
pub const ln_ = ops.ln_;
// pub const Log = numops.Log;
// pub const log = numops.log;
// pub const log_ = numops.log_;

// Power functions
pub const Pow = ops.Pow;
pub const pow = ops.pow;
pub const pow_ = ops.pow_;
pub const Sqrt = ops.Sqrt;
pub const sqrt = ops.sqrt;
pub const sqrt_ = ops.sqrt_;
pub const Cbrt = ops.Cbrt;
pub const cbrt = ops.cbrt;
pub const cbrt_ = ops.cbrt_;
pub const Hypot = ops.Hypot;
pub const hypot = ops.hypot;
pub const hypot_ = ops.hypot_;

// Trigonometric functions
pub const Sin = ops.Sin;
pub const sin = ops.sin;
pub const sin_ = ops.sin_;
pub const Cos = ops.Cos;
pub const cos = ops.cos;
pub const cos_ = ops.cos_;
pub const Tan = ops.Tan;
pub const tan = ops.tan;
pub const tan_ = ops.tan_;
pub const Asin = ops.Asin;
pub const asin = ops.asin;
pub const asin_ = ops.asin_;
pub const Acos = ops.Acos;
pub const acos = ops.acos;
pub const acos_ = ops.acos_;
pub const Atan = ops.Atan;
pub const atan = ops.atan;
pub const atan_ = ops.atan_;
pub const Atan2 = ops.Atan2;
pub const atan2 = ops.atan2;
pub const atan2_ = ops.atan2_;

// Hyperbolic functions
pub const Sinh = ops.Sinh;
pub const sinh = ops.sinh;
pub const sinh_ = ops.sinh_;
pub const Cosh = ops.Cosh;
pub const cosh = ops.cosh;
pub const cosh_ = ops.cosh_;
pub const Tanh = ops.Tanh;
pub const tanh = ops.tanh;
pub const tanh_ = ops.tanh_;
pub const Asinh = ops.Asinh;
pub const asinh = ops.asinh;
pub const asinh_ = ops.asinh_;
pub const Acosh = ops.Acosh;
pub const acosh = ops.acosh;
pub const acosh_ = ops.acosh_;
pub const Atanh = ops.Atanh;
pub const atanh = ops.atanh;
pub const atanh_ = ops.atanh_;

// Special functions
pub const Erf = ops.Erf;
pub const erf = ops.erf;
pub const erf_ = ops.erf_;
pub const Erfc = ops.Erfc;
pub const erfc = ops.erfc;
pub const erfc_ = ops.erfc_;
pub const Gamma = ops.Gamma;
pub const gamma = ops.gamma;
pub const gamma_ = ops.gamma_;
pub const Lgamma = ops.Lgamma;
pub const lgamma = ops.lgamma;
pub const lgamma_ = ops.lgamma_;
