//! Namespace for float operations.

const types = @import("types.zig");
const Cmp = types.Cmp;

// Constant functions
pub inline fn pi(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .float)
        @compileError("zml.float.pi: T must be a float type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return 3.1415926535897932384626433832795028841971;
}

pub inline fn pi_2(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .float)
        @compileError("zml.float.pi_2: T must be a float type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return 1.5707963267948966192313216916397514420985;
}

pub inline fn pi_4(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .float)
        @compileError("zml.float.pi_4: T must be a float type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return 0.7853981633974483096156608458198757212951;
}

pub inline fn e(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .float)
        @compileError("zml.float.e: T must be a float type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return 2.7182818284590452353602874713526624977572;
}

pub inline fn ln2(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .float)
        @compileError("zml.float.ln2: T must be a float type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return 0.6931471805599453094172321214581765680755;
}

pub inline fn log10e(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .float)
        @compileError("zml.float.log10e: T must be a float type, got \n\tT: " ++ @typeName(T) ++ "\n");

    return 0.4342944819032518276511289189166050822944;
}

// Basic functions
pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.add: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn add(
    x: anytype,
    y: anytype,
) Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = Add(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x) + types.scast(R, y);
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.sub: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn sub(
    x: anytype,
    y: anytype,
) Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = Sub(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x) - types.scast(R, y);
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.mul: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn mul(
    x: anytype,
    y: anytype,
) Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = Mul(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x) * types.scast(R, y);
}

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or !types.isNumeric(Z) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or !types.numericType(Z).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float and types.numericType(Z) != .float))
        @compileError("zml.float.fma: at least one of x, y or z must be a float, the others must be bool, int or float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\tz: " ++ @typeName(Z) ++ "\n");

    return types.Coerce(X, types.Coerce(Y, Z));
}

pub inline fn fma(
    x: anytype,
    y: anytype,
    z: anytype,
) Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const R: type = Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z));

    return @mulAdd(R, types.scast(R, x), types.scast(R, y), types.scast(R, z));
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.div: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn div(
    x: anytype,
    y: anytype,
) Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = Div(@TypeOf(x), @TypeOf(y));

    return types.scast(R, x) / types.scast(R, y);
}

pub inline fn cmp(
    x: anytype,
    y: anytype,
) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.cmp: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    if (types.scast(C, x) < types.scast(C, y)) return .lt;
    if (types.scast(C, x) > types.scast(C, y)) return .gt;
    return .eq;
}

pub inline fn eq(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.eq: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) == types.scast(C, y);
}

pub inline fn ne(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.ne: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) != types.scast(C, y);
}

pub inline fn lt(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.lt: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) < types.scast(C, y);
}

pub inline fn le(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.le: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) <= types.scast(C, y);
}

pub inline fn gt(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.gt: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) > types.scast(C, y);
}

pub inline fn ge(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.ge: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) >= types.scast(C, y);
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.max: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn max(
    x: anytype,
    y: anytype,
) Max(@TypeOf(x), @TypeOf(y)) {
    const R: type = Max(@TypeOf(x), @TypeOf(y));

    return if (types.scast(R, x) > types.scast(R, y)) types.scast(R, x) else types.scast(R, y);
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.float) or !types.numericType(Y).le(.float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("zml.float.min: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn min(
    x: anytype,
    y: anytype,
) Min(@TypeOf(x), @TypeOf(y)) {
    const R: type = Min(@TypeOf(x), @TypeOf(y));

    return if (types.scast(R, x) < types.scast(R, y)) types.scast(R, x) else types.scast(R, y);
}

// Basic operations
pub const abs = @import("float/abs.zig").abs;
// pub const fmod = @import("float/fmod.zig").fmod; // to implement
// pub const remainder = @import("float/remainder.zig").remainder; // to implement
// pub const remquo = @import("float/remquo.zig").remquo; // to implement
// pub const fdim = @import("float/fdim.zig").fdim; // to implement

// Exponential functions
pub const Exp = @import("float/exp.zig").Exp;
pub const exp = @import("float/exp.zig").exp;
pub const Exp2 = @import("float/exp2.zig").Exp2;
pub const exp2 = @import("float/exp2.zig").exp2;
pub const Expm1 = @import("float/expm1.zig").Expm1;
pub const expm1 = @import("float/expm1.zig").expm1;
pub const Log = @import("float/log.zig").Log;
pub const log = @import("float/log.zig").log;
pub const Log10 = @import("float/log10.zig").Log10;
pub const log10 = @import("float/log10.zig").log10;
pub const Log2 = @import("float/log2.zig").Log2;
pub const log2 = @import("float/log2.zig").log2;
pub const Log1p = @import("float/log1p.zig").Log1p;
pub const log1p = @import("float/log1p.zig").log1p;

// Power functions
pub const Pow = @import("float/pow.zig").Pow;
pub const pow = @import("float/pow.zig").pow;
pub const Sqrt = @import("float/sqrt.zig").Sqrt;
pub const sqrt = @import("float/sqrt.zig").sqrt;
pub const Cbrt = @import("float/cbrt.zig").Cbrt;
pub const cbrt = @import("float/cbrt.zig").cbrt;
pub const Hypot = @import("float/hypot.zig").Hypot;
pub const hypot = @import("float/hypot.zig").hypot;

// Trigonometric functions
pub const Sin = @import("float/sin.zig").Sin;
pub const sin = @import("float/sin.zig").sin;
pub const Cos = @import("float/cos.zig").Cos;
pub const cos = @import("float/cos.zig").cos;
pub const Tan = @import("float/tan.zig").Tan;
pub const tan = @import("float/tan.zig").tan;
pub const Asin = @import("float/asin.zig").Asin;
pub const asin = @import("float/asin.zig").asin;
pub const Acos = @import("float/acos.zig").Acos;
pub const acos = @import("float/acos.zig").acos;
pub const Atan = @import("float/atan.zig").Atan;
pub const atan = @import("float/atan.zig").atan;
pub const Atan2 = @import("float/atan2.zig").Atan2;
pub const atan2 = @import("float/atan2.zig").atan2;
pub const Sincos = @import("float/sincos.zig").Sincos;
pub const sincos = @import("float/sincos.zig").sincos;

// Hyperbolic functions
pub const Sinh = @import("float/sinh.zig").Sinh;
pub const sinh = @import("float/sinh.zig").sinh;
pub const Cosh = @import("float/cosh.zig").Cosh;
pub const cosh = @import("float/cosh.zig").cosh;
pub const Tanh = @import("float/tanh.zig").Tanh;
pub const tanh = @import("float/tanh.zig").tanh;
pub const Asinh = @import("float/asinh.zig").Asinh;
pub const asinh = @import("float/asinh.zig").asinh;
pub const Acosh = @import("float/acosh.zig").Acosh;
pub const acosh = @import("float/acosh.zig").acosh;
pub const Atanh = @import("float/atanh.zig").Atanh;
pub const atanh = @import("float/atanh.zig").atanh;

// Error and gamma functions
pub const Erf = @import("float/erf.zig").Erf;
pub const erf = @import("float/erf.zig").erf;
pub const Erfc = @import("float/erfc.zig").Erfc;
pub const erfc = @import("float/erfc.zig").erfc;
pub const Gamma = @import("float/gamma.zig").Gamma;
pub const gamma = @import("float/gamma.zig").gamma;
pub const Lgamma = @import("float/lgamma.zig").Lgamma;
pub const lgamma = @import("float/lgamma.zig").lgamma;

// Bessel functions
// pub const j0 = @import("float/j0.zig").j0; // to implement
// pub const j1 = @import("float/j1.zig").j1; // to implement
// pub const jn = @import("float/jn.zig").jn; // to implement
// pub const y0 = @import("float/y0.zig").y0; // to implement
// pub const y1 = @import("float/y1.zig").y1; // to implement
// pub const yn = @import("float/yn.zig").yn; // to implement

// Nearest integer floating-point operations
pub const ceil = @import("float/ceil.zig").ceil;
pub const floor = @import("float/floor.zig").floor;
pub const trunc = @import("float/trunc.zig").trunc;
pub const round = @import("float/round.zig").round;
// pub const nearbyint = @import("float/nearbyint.zig").nearbyint; // to implement
pub const rint = @import("float/rint.zig").rint;

// Floating-point manipulation functions
pub const Frexp = @import("float/frexp.zig").Frexp;
pub const frexp = @import("float/frexp.zig").frexp;
pub const Ldexp = @import("float/ldexp.zig").Ldexp;
pub const ldexp = @import("float/ldexp.zig").ldexp;
// pub const modf = @import("float/modf.zig").modf; // to implement
pub const Scalbn = @import("float/scalbn.zig").Scalbn;
pub const scalbn = @import("float/scalbn.zig").scalbn;
// pub const ilogb = @import("float/ilogb.zig").ilogb; // to implement
// pub const logb = @import("float/logb.zig").logb; // to implement
// pub const nextafter = @import("float/nextafter.zig").nextafter; // to implement
// pub const nexttoward = @import("float/nexttoward.zig").nexttoward; // to implement
pub const copysign = @import("float/copysign.zig").copysign;

pub const Error = error{
    NotFinite,
};
