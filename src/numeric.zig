//! Namespace for numeric types and operations.

const constants = @import("numeric/constants.zig");
const ops = @import("numeric/ops.zig");

// Constants
pub const highest = constants.highest;
pub const lowest = constants.lowest;
pub const smallest = constants.smallest;
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
pub const absInto = ops.absInto;
pub const Abs1 = ops.Abs1;
pub const abs1 = ops.abs1;
pub const abs1Into = ops.abs1Into;
pub const Abs2 = ops.Abs2;
pub const abs2 = ops.abs2;
pub const abs2Into = ops.abs2Into;
pub const Neg = ops.Neg;
pub const neg = ops.neg;
pub const negInto = ops.negInto;
pub const Re = ops.Re;
pub const re = ops.re;
pub const Im = ops.Im;
pub const im = ops.im;
pub const Conj = ops.Conj;
pub const conj = ops.conj;
pub const conjInto = ops.conjInto;
pub const Sign = ops.Sign;
pub const sign = ops.sign;
// pub const copysign = numops.copysign;

// Arithmetic operations
pub const Add = ops.Add;
pub const add = ops.add;
pub const addInto = ops.addInto;
pub const Sub = ops.Sub;
pub const sub = ops.sub;
pub const subInto = ops.subInto;
pub const Mul = ops.Mul;
pub const mul = ops.mul;
pub const mulInto = ops.mulInto;
pub const Fma = ops.Fma;
pub const fma = ops.fma;
pub const fmaInto = ops.fmaInto;
pub const Div = ops.Div;
pub const div = ops.div;
pub const divInto = ops.divInto;

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
pub const maxInto = ops.maxInto;
pub const Min = ops.Min;
pub const min = ops.min;
pub const minInto = ops.minInto;

// Exponential functions
pub const Exp = ops.Exp;
pub const exp = ops.exp;
pub const expInto = ops.expInto;
pub const Ln = ops.Ln;
pub const ln = ops.ln;
pub const lnInto = ops.lnInto;
// pub const Log = numops.Log;
// pub const log = numops.log;
// pub const logInto = numops.logInto;

// Power functions
pub const Pow = ops.Pow;
pub const pow = ops.pow;
pub const powInto = ops.powInto;
pub const Sqrt = ops.Sqrt;
pub const sqrt = ops.sqrt;
pub const sqrtInto = ops.sqrtInto;
pub const Cbrt = ops.Cbrt;
pub const cbrt = ops.cbrt;
pub const cbrtInto = ops.cbrtInto;
pub const Hypot = ops.Hypot;
pub const hypot = ops.hypot;
pub const hypotInto = ops.hypotInto;

// Trigonometric functions
pub const Sin = ops.Sin;
pub const sin = ops.sin;
pub const sinInto = ops.sinInto;
pub const Cos = ops.Cos;
pub const cos = ops.cos;
pub const cosInto = ops.cosInto;
pub const Tan = ops.Tan;
pub const tan = ops.tan;
pub const tanInto = ops.tanInto;
pub const Asin = ops.Asin;
pub const asin = ops.asin;
pub const asinInto = ops.asinInto;
pub const Acos = ops.Acos;
pub const acos = ops.acos;
pub const acosInto = ops.acosInto;
pub const Atan = ops.Atan;
pub const atan = ops.atan;
pub const atanInto = ops.atanInto;
pub const Atan2 = ops.Atan2;
pub const atan2 = ops.atan2;
pub const atan2Into = ops.atan2Into;

// Hyperbolic functions
pub const Sinh = ops.Sinh;
pub const sinh = ops.sinh;
pub const sinhInto = ops.sinhInto;
pub const Cosh = ops.Cosh;
pub const cosh = ops.cosh;
pub const coshInto = ops.coshInto;
pub const Tanh = ops.Tanh;
pub const tanh = ops.tanh;
pub const tanhInto = ops.tanhInto;
pub const Asinh = ops.Asinh;
pub const asinh = ops.asinh;
pub const asinhInto = ops.asinhInto;
pub const Acosh = ops.Acosh;
pub const acosh = ops.acosh;
pub const acoshInto = ops.acoshInto;
pub const Atanh = ops.Atanh;
pub const atanh = ops.atanh;
pub const atanhInto = ops.atanhInto;

// Special functions
pub const Erf = ops.Erf;
pub const erf = ops.erf;
pub const erfInto = ops.erfInto;
pub const Erfc = ops.Erfc;
pub const erfc = ops.erfc;
pub const erfcInto = ops.erfcInto;
pub const Gamma = ops.Gamma;
pub const gamma = ops.gamma;
pub const gammaInto = ops.gammaInto;
pub const Lgamma = ops.Lgamma;
pub const lgamma = ops.lgamma;
pub const lgammaInto = ops.lgammaInto;
