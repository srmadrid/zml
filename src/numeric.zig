//! Namespace for numeric types and operations.

const int = @import("int.zig");
const float = @import("float.zig");
const dyadic = @import("dyadic.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const complex = @import("complex.zig");

const numops = @import("numeric/ops.zig");

// Basic operations
pub const Abs = numops.Abs;
pub const abs = numops.abs;
pub const abs_ = numops.abs_;
pub const Abs1 = numops.Abs1;
pub const abs1 = numops.abs1;
pub const abs1_ = numops.abs1_;
pub const Abs2 = numops.Abs2;
pub const abs2 = numops.abs2;
pub const abs2_ = numops.abs2_;
pub const Neg = numops.Neg;
pub const neg = numops.neg;
pub const neg_ = numops.neg_;
pub const Re = numops.Re;
pub const re = numops.re;
pub const Im = numops.Im;
pub const im = numops.im;
pub const Conj = numops.Conj;
pub const conj = numops.conj;
pub const conj_ = numops.conj_;
pub const copysign = numops.copysign;

// Arithmetic operations
pub const Add = numops.Add;
pub const add = numops.add;
pub const add_ = numops.add_;
pub const Sub = numops.Sub;
pub const sub = numops.sub;
pub const sub_ = numops.sub_;
pub const Mul = numops.Mul;
pub const mul = numops.mul;
pub const mul_ = numops.mul_;
pub const Div = numops.Div;
pub const div = numops.div;
pub const div_ = numops.div_;

// Comparison operations
pub const cmp = numops.cmp;
pub const Eq = numops.Eq;
pub const eq = numops.eq;
pub const eq_ = numops.eq_;
pub const Ne = numops.Ne;
pub const ne = numops.ne;
pub const ne_ = numops.ne_;
pub const Lt = numops.Lt;
pub const lt = numops.lt;
pub const lt_ = numops.lt_;
pub const Le = numops.Le;
pub const le = numops.le;
pub const le_ = numops.le_;
pub const Gt = numops.Gt;
pub const gt = numops.gt;
pub const gt_ = numops.gt_;
pub const Ge = numops.Ge;
pub const ge = numops.ge;
pub const ge_ = numops.ge_;
pub const Max = numops.Max;
pub const max = numops.max;
pub const max_ = numops.max_;
pub const Min = numops.Min;
pub const min = numops.min;
pub const min_ = numops.min_;

// Exponential functions
pub const Exp = numops.Exp;
pub const exp = numops.exp;
pub const exp_ = numops.exp_;
pub const Exp2 = numops.Exp2;
pub const exp2 = numops.exp2;
pub const exp2_ = numops.exp2_;
pub const Log = numops.Log;
pub const log = numops.log;
pub const log_ = numops.log_;
pub const Log10 = numops.Log10;
pub const log10 = numops.log10;
pub const log10_ = numops.log10_;
pub const Log2 = numops.Log2;
pub const log2 = numops.log2;
pub const log2_ = numops.log2_;

// Power functions
pub const Pow = numops.Pow;
pub const pow = numops.pow;
pub const pow_ = numops.pow_;
pub const Sqrt = numops.Sqrt;
pub const sqrt = numops.sqrt;
pub const sqrt_ = numops.sqrt_;
pub const Cbrt = numops.Cbrt;
pub const cbrt = numops.cbrt;
pub const cbrt_ = numops.cbrt_;
pub const Hypot = numops.Hypot;
pub const hypot = numops.hypot;
pub const hypot_ = numops.hypot_;

// Trigonometric functions
pub const Sin = numops.Sin;
pub const sin = numops.sin;
pub const sin_ = numops.sin_;
pub const Cos = numops.Cos;
pub const cos = numops.cos;
pub const cos_ = numops.cos_;
pub const Tan = numops.Tan;
pub const tan = numops.tan;
pub const tan_ = numops.tan_;
pub const Asin = numops.Asin;
pub const asin = numops.asin;
pub const asin_ = numops.asin_;
pub const Acos = numops.Acos;
pub const acos = numops.acos;
pub const acos_ = numops.acos_;
pub const Atan = numops.Atan;
pub const atan = numops.atan;
pub const atan_ = numops.atan_;
pub const Atan2 = numops.Atan2;
pub const atan2 = numops.atan2;
pub const atan2_ = numops.atan2_;

// Hyperbolic functions
pub const Sinh = numops.Sinh;
pub const sinh = numops.sinh;
pub const sinh_ = numops.sinh_;
pub const Cosh = numops.Cosh;
pub const cosh = numops.cosh;
pub const cosh_ = numops.cosh_;
pub const Tanh = numops.Tanh;
pub const tanh = numops.tanh;
pub const tanh_ = numops.tanh_;
pub const Asinh = numops.Asinh;
pub const asinh = numops.asinh;
pub const asinh_ = numops.asinh_;
pub const Acosh = numops.Acosh;
pub const acosh = numops.acosh;
pub const acosh_ = numops.acosh_;
pub const Atanh = numops.Atanh;
pub const atanh = numops.atanh;
pub const atanh_ = numops.atanh_;

// Special functions
pub const Erf = numops.Erf;
pub const erf = numops.erf;
pub const erf_ = numops.erf_;
pub const Erfc = numops.Erfc;
pub const erfc = numops.erfc;
pub const erfc_ = numops.erfc_;
pub const Gamma = numops.Gamma;
pub const gamma = numops.gamma;
pub const gamma_ = numops.gamma_;
pub const Lgamma = numops.Lgamma;
pub const lgamma = numops.lgamma;
pub const lgamma_ = numops.lgamma_;
