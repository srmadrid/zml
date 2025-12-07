//! Namespace for float operations.

const types = @import("types.zig");
const scast = types.scast;
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const Order = types.Order;

// Constant functions
pub inline fn pi(comptime T: type) T {
    comptime if (types.numericType(T) != .float)
        @compileError("T must be a float type");

    return 3.1415926535897932384626433832795028841971;
}

pub inline fn pi_2(comptime T: type) T {
    comptime if (types.numericType(T) != .float)
        @compileError("T must be a float type");

    return 1.5707963267948966192313216916397514420985;
}

pub inline fn pi_4(comptime T: type) T {
    comptime if (types.numericType(T) != .float)
        @compileError("T must be a float type");

    return 0.7853981633974483096156608458198757212951;
}

pub inline fn e(comptime T: type) T {
    comptime if (types.numericType(T) != .float)
        @compileError("T must be a float type");

    return 2.7182818284590452353602874713526624977572;
}

pub inline fn ln2(comptime T: type) T {
    comptime if (types.numericType(T) != .float)
        @compileError("T must be a float type");

    return 0.6931471805599453094172321214581765680755;
}

pub inline fn log10e(comptime T: type) T {
    comptime if (types.numericType(T) != .float)
        @compileError("T must be a float type");

    return 0.4342944819032518276511289189166050822944;
}

// Basic functions
pub inline fn add(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("float.add requires at least one of x or y to be a float, the other must be a bool, int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return scast(C, x) + scast(C, y);
}

pub inline fn sub(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("float.sub requires at least one of x or y to be a float, the other must be a bool, int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return scast(C, x) - scast(C, y);
}

pub inline fn mul(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("float.mul requires at least one of x or y to be a float, the other must be a bool, int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return scast(C, x) * scast(C, y);
}

pub inline fn div(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float) or
        (types.numericType(X) != .float and types.numericType(Y) != .float))
        @compileError("float.div requires at least one of x or y to be a float, the other must be a bool, int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return scast(C, x) / scast(C, y);
}

pub inline fn cmp(
    x: anytype,
    y: anytype,
) Order {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.cmp requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (x < y) return .lt;
    if (x > y) return .gt;
    return .eq;
}

pub inline fn eq(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.eq requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return types.scast(Coerce(X, Y), x) == types.scast(Coerce(X, Y), y);
}

pub inline fn ne(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.neq requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return types.scast(Coerce(X, Y), x) != types.scast(Coerce(X, Y), y);
}

pub inline fn lt(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.lt requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return types.scast(Coerce(X, Y), x) < types.scast(Coerce(X, Y), y);
}

pub inline fn le(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.le requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return types.scast(Coerce(X, Y), x) <= types.scast(Coerce(X, Y), y);
}

pub inline fn gt(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.gt requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return types.scast(Coerce(X, Y), x) > types.scast(Coerce(X, Y), y);
}

pub inline fn ge(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.ge requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return types.scast(Coerce(X, Y), x) >= types.scast(Coerce(X, Y), y);
}

pub inline fn max(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.ge requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return if (x > y) scast(C, x) else scast(C, y);
}

pub inline fn min(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (!(types.numericType(X) == .int and types.numericType(X) == .float) and
        (types.numericType(Y) == .int and types.numericType(Y) == .float) and
        (types.numericType(X) == .float and types.numericType(Y) == .float))
        @compileError("float.min requires at least one of x or y to be a float, the other must be an int or float, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return if (x < y) scast(C, x) else scast(C, y);
}

// Basic operations
pub const abs = @import("float/abs.zig").abs;
// pub const fmod = @import("float/fmod.zig").fmod; // to implement
// pub const remainder = @import("float/remainder.zig").remainder; // to implement
// pub const remquo = @import("float/remquo.zig").remquo; // to implement
// pub const fdim = @import("float/fdim.zig").fdim; // to implement

// Exponential functions
pub const exp = @import("float/exp.zig").exp;
pub const exp10 = @import("float/exp10.zig").exp10;
pub const exp2 = @import("float/exp2.zig").exp2;
pub const exp10m1 = @import("float/exp10m1.zig").exp10m1;
pub const exp2m1 = @import("float/exp2m1.zig").exp2m1;
pub const expm1 = @import("float/expm1.zig").expm1;
pub const log = @import("float/log.zig").log;
pub const log10 = @import("float/log10.zig").log10;
pub const log2 = @import("float/log2.zig").log2;
pub const log10p1 = @import("float/log10p1.zig").log10p1;
pub const log2p1 = @import("float/log2p1.zig").log2p1;
pub const log1p = @import("float/log1p.zig").log1p;

// Power functions
pub const pow = @import("float/pow.zig").pow;
pub const sqrt = @import("float/sqrt.zig").sqrt;
pub const cbrt = @import("float/cbrt.zig").cbrt;
pub const hypot = @import("float/hypot.zig").hypot;

// Trigonometric functions
pub const sin = @import("float/sin.zig").sin;
pub const cos = @import("float/cos.zig").cos;
pub const tan = @import("float/tan.zig").tan;
pub const asin = @import("float/asin.zig").asin;
pub const acos = @import("float/acos.zig").acos;
pub const atan = @import("float/atan.zig").atan;
pub const atan2 = @import("float/atan2.zig").atan2;
pub const sincos = @import("float/sincos.zig").sincos;

// Hyperbolic functions
pub const sinh = @import("float/sinh.zig").sinh;
pub const cosh = @import("float/cosh.zig").cosh;
pub const tanh = @import("float/tanh.zig").tanh;
pub const asinh = @import("float/asinh.zig").asinh;
pub const acosh = @import("float/acosh.zig").acosh;
pub const atanh = @import("float/atanh.zig").atanh;

// Error and gamma functions
pub const erf = @import("float/erf.zig").erf;
pub const erfc = @import("float/erfc.zig").erfc;
pub const gamma = @import("float/gamma.zig").gamma;
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
// pub const lrint = @import("float/lrint.zig").lrint; // to implement
// pub const llrint = @import("float/llrint.zig").llrint; // to implement. Join the rint's and choose one depending on input integer type?

// Floating-point manipulation functions
pub const frexp = @import("float/frexp.zig").frexp;
pub const ldexp = @import("float/ldexp.zig").ldexp;
// pub const modf = @import("float/modf.zig").modf; // to implement
pub const scalbn = @import("float/scalbn.zig").scalbn;
// pub const ilogb = @import("float/ilogb.zig").ilogb; // to implement
// pub const logb = @import("float/logb.zig").logb; // to implement
// pub const nextafter = @import("float/nextafter.zig").nextafter; // to implement
// pub const nexttoward = @import("float/nexttoward.zig").nexttoward; // to implement
pub const copysign = @import("float/copysign.zig").copysign;
