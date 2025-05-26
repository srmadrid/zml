const types = @import("types.zig");
const scast = types.scast;
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;

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
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.add requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    return scast(C, left) + scast(C, right);
}

pub inline fn add_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("float.add_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.add_ requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    out.* = scast(O, scast(C, left) + scast(C, right));
}

pub inline fn sub(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.sub requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    return scast(C, left) - scast(C, right);
}

pub inline fn sub_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("float.sub_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.sub_ requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    out.* = scast(O, scast(C, left) - scast(C, right));
}

pub inline fn mul(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.mul requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    return scast(C, left) * scast(C, right);
}

pub inline fn mul_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("float.mul_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.mul_ requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    out.* = scast(O, scast(C, left) * scast(C, right));
}

pub inline fn div(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.div requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    return scast(C, left) / scast(C, right);
}

pub inline fn div_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("float.div_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float) or
        (types.numericType(L) != .float and types.numericType(R) != .float))
        @compileError("float.div_ requires at least one of left or right to be a float, the other must be a bool, int or float, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    out.* = scast(O, scast(C, left) / scast(C, right));
}

pub const abs = @import("float/abs.zig").abs;
// pub const fmod = @import("float/fmod.zig").fmod; // to implement
// pub const remainder = @import("float/remainder.zig").remainder; // to implement
// pub const remquo = @import("float/remquo.zig").remquo; // to implement
// pub const fdim = @import("float/fdim.zig").fdim; // to implement

// Exponential functions
pub const exp = @import("float/exp.zig").exp; // to implement: f80. 5/768 tests fail: 1 for f32, 1 for f64, 3 for f128
pub const exp10 = @import("float/exp10.zig").exp10; // to implement: f80. 54/650 tests fail: 54 for f128
pub const exp2 = @import("float/exp2.zig").exp2; // to implement: f80. 5/670 tests fail: 1 for f64, 4 for f128
pub const exp10m1 = @import("float/exp10m1.zig").exp10m1; // 235/977 tests fail: 28 for f32, 68 for f64, 55 for f80, 84 for f128
pub const exp2m1 = @import("float/exp2m1.zig").exp2m1; // 189/957 tests fail: 24 for f32, 49 for f64, 48 for f80, 68 for f128
pub const expm1 = @import("float/expm1.zig").expm1; // to implement: f80. 28/507 tests fail: 12 for f64, 16 for f128
pub const log = @import("float/log.zig").log; // to implement: f80. 8/256 tests fail: 8 for f128
pub const log10 = @import("float/log10.zig").log10; // to implement: f80. 33/264 tests fail: 16 for f64, 17 for f128
pub const log2 = @import("float/log2.zig").log2; // to implement: f80. 35/326 tests fail: 5 for f32, 30 for f128
pub const log10p1 = @import("float/log10p1.zig").log10p1; // 240/749 tests fail: 38 for f32, 55 for f64, 59 for f80, 88 for f128
pub const log2p1 = @import("float/log2p1.zig").log2p1; // 122/473 tests fail: 16 for f32, 15 for f64, 25 for f80, 66 for f128
pub const log1p = @import("float/log1p.zig").log1p; // to implement: f80. 30/421 tests fail: 8 for f64, 22 for f128

// Power functions
pub const pow = @import("float/pow.zig").pow; // to implement: f80. 41/6896 tests fail: 2 for f32, 2 for f64, 2 for f80, 35 for f128
pub const sqrt = @import("float/sqrt.zig").sqrt; // 58/655 tests fail: 1 for f80, 57 for f128
pub const cbrt = @import("float/cbrt.zig").cbrt; // to implement: f80. 23/233 tests fail: 21 for f64, 2 for f128
pub const hypot = @import("float/hypot.zig").hypot; // to implement: f80. 2/2281 tests fail: 1 for f32, 1 for f128

// Trigonometric functions
pub const sin = @import("float/sin.zig").sin; // to implement: f80. 20/609 tests fail: 1 for f32, 6 for f64, 13 for f128
pub const cos = @import("float/cos.zig").cos; // to implement: f80. 16/532 tests fail: 2 for f32, 2 for f64, 12 for f128
pub const tan = @import("float/tan.zig").tan; // to implement: f80. 3/546 tests fail: 3 for f128
pub const asin = @import("float/asin.zig").asin; // to implement: f80. 3/385 tests fail: 1 for f64, 2 for f128
pub const acos = @import("float/acos.zig").acos; // to implement: f80. 2/494 tests fail: 1 for f64, 1 for f128
pub const atan = @import("float/atan.zig").atan; // to implement: f80. 3/228 tests fail: 1 for f64, 2 for f128
pub const atan2 = @import("float/atan2.zig").atan2; // to implement: f80. 107/2672 tests fail: 2 for f80, 105 for f128
pub const sincos = @import("float/sincos.zig").sincos; // to implement: f80. 8/674 tests fail: 1 for f64, 7 for f128
pub const sinpi = @import("float/sinpi.zig").sinpi; // 154/1572 tests fail: 28 for f32, 25 for f64, 38 for f80, 63 for f128
pub const cospi = @import("float/cospi.zig").cospi; // 266/1468 tests fail: 43 for f32, 39 for f64, 70 for f80, 114 for f128
pub const tanpi = @import("float/tanpi.zig").tanpi; // 535/1436 tests fail: 58 for f32, 131 for f64, 128 for f80, 218 for f128
pub const asinpi = @import("float/asinpi.zig").asinpi; // 121/437 tests fail: 25 for f32, 20 f64, 30 for f80, 46 for f128
pub const acospi = @import("float/acospi.zig").acospi; // 131/546 tests fail: 17 for f32, 29 for f64, 29 for f80, 56 for f128
pub const atanpi = @import("float/atanpi.zig").atanpi; // 51/272 tests fail: 8 for f32, 11 for f64, 9 for f80, 23 for f128
pub const atan2pi = @import("float/atan2pi.zig").atan2pi; // 177/2413: 7 for f32, 30 for f64, 45 for f80, 95 for f128

// Hyperbolic functions
pub const sinh = @import("float/sinh.zig").sinh; // to implement: f80. 67/577 tests fail: 29 for f64, 38 for f128
pub const cosh = @import("float/cosh.zig").cosh; // to implement: f80. 52/577 tests fail: 12 for f64, 40 for f128
pub const tanh = @import("float/tanh.zig").tanh; // to implement: f80. 40/452 tests fail: 15 for f64, 25 for f128
pub const asinh = @import("float/asinh.zig").asinh; // to implement: f80. 51/507 tests fail: 16 for f64, 35 for f128
pub const acosh = @import("float/acosh.zig").acosh; // to implement: f80. 51/349 tests fail: 23 for f64, 28 for f128
pub const atanh = @import("float/atanh.zig").atanh; // to implement: f80. 48/566 tests fail: 17 for f64, 31 for f128

// Error and gamma functions
pub const erf = @import("float/erf.zig").erf; // to implement: f80. 11/532 tests fail: 1 for f64, 10 for f128
pub const erfc = @import("float/erfc.zig").erfc; // to implement: f80. 142/705 tests fail: 36 for f64, 22 for f80, 84 for f128
pub const gamma = @import("float/gamma.zig").gamma; // to implement: f80. 943/4212 tests fail: 308 for f64, 635 for f128
pub const lgamma = @import("float/lgamma.zig").lgamma; // to implement: f80. 874/3800 tests fail: 303 for f64, 571 for f128

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
