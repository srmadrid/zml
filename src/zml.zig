pub const types = @import("types.zig");
pub const scast = types.scast;
pub const cast = types.cast;
pub const Cmp = types.Cmp;
pub const Layout = types.Layout;
pub const Uplo = types.Uplo;
pub const Diag = types.Diag;
pub const IterationOrder = types.IterationOrder;

pub const int = @import("int.zig");
pub const float = @import("float.zig");
pub const dyadic = @import("dyadic.zig");
pub const Dyadic = dyadic.Dyadic;
pub const cfloat = @import("cfloat.zig");
pub const cf16 = cfloat.cf16;
pub const cf32 = cfloat.cf32;
pub const cf64 = cfloat.cf64;
pub const cf80 = cfloat.cf80;
pub const cf128 = cfloat.cf128;
pub const comptime_cfloat = cfloat.comptime_cfloat;
pub const integer = @import("integer.zig");
pub const Integer = integer.Integer;
pub const rational = @import("rational.zig");
pub const Rational = rational.Rational;
pub const real = @import("real.zig");
pub const Real = real.Real;
pub const complex = @import("complex.zig");
pub const Complex = complex.Complex;

const constants = @import("constants.zig");
pub const zero = constants.zero;
pub const one = constants.one;
pub const two = constants.two;

const ops = @import("ops.zig");

// Arithmetic operations
pub const Add = ops.Add;
pub const add = ops.add;
pub const add_ = ops.add_;
pub const Sub = ops.Sub;
pub const sub = ops.sub;
pub const sub_ = ops.sub_;
pub const Mul = ops.Mul;
pub const mul = ops.mul;
pub const mul_ = ops.mul_;
pub const Div = ops.Div;
pub const div = ops.div;
pub const div_ = ops.div_;

// Comparison operations
pub const Eq = ops.Eq;
pub const eq = ops.eq;
pub const eq_ = ops.eq_;
pub const Ne = ops.Ne;
pub const ne = ops.ne;
pub const ne_ = ops.ne_;
pub const Lt = ops.Lt;
pub const lt = ops.lt;
pub const lt_ = ops.lt_;
pub const Le = ops.Le;
pub const le = ops.le;
pub const le_ = ops.le_;
pub const Gt = ops.Gt;
pub const gt = ops.gt;
pub const gt_ = ops.gt_;
pub const Ge = ops.Ge;
pub const ge = ops.ge;
pub const ge_ = ops.ge_;
pub const Max = ops.Max;
pub const max = ops.max;
pub const max_ = ops.max_;
pub const Min = ops.Min;
pub const min = ops.min;
pub const min_ = ops.min_;

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

// Exponential functions
pub const Exp = ops.Exp;
pub const exp = ops.exp;
pub const exp_ = ops.exp_;
pub const Exp2 = ops.Exp2;
pub const exp2 = ops.exp2;
pub const exp2_ = ops.exp2_;
pub const Log = ops.Log;
pub const log = ops.log;
pub const log_ = ops.log_;
pub const Log2 = ops.Log2;
pub const log2 = ops.log2;
pub const log2_ = ops.log2_;
pub const Log10 = ops.Log10;
pub const log10 = ops.log10;
pub const log10_ = ops.log10_;

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

// Error and gamma functions
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

// Nearest integer operations
pub const ceil = ops.ceil;
pub const ceil_ = ops.ceil_;

pub const set = ops.set;
pub const deinit = ops.deinit;

// Domain namespaces
pub const numeric = @import("numeric.zig");
pub const vector = @import("vector.zig");
pub const matrix = @import("matrix.zig");
pub const array = @import("array.zig");

pub const linalg = @import("linalg.zig");
pub const autodiff = @import("autodiff.zig");

// Symbolic system.
//pub const Expression = @import("expression/expression.zig").Expression;
//pub const Symbol = @import("symbol.zig").Symbol;
//pub const Element = @import("element.zig").Element;
//pub const Variable = @import("variable.zig").Variable;
//pub const Set = @import("set.zig").Set;
