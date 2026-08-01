//! Namespace for float operations.

const float = @This();

const std = @import("std");

const meta = @import("meta.zig");
const numeric = @import("numeric.zig");

pub const Coerce = @import("float/coerce.zig").Coerce;

// Constant functions
pub fn pi(comptime Float: type) Float {
    comptime if (!meta.isNumeric(Float) or meta.numericType(Float) != .float)
        @compileError("zsl.float.pi: Float must be a float type, got \n\nFloat = " ++ @typeName(Float) ++ "\n");

    return 3.1415926535897932384626433832795028841971;
}

pub fn e(comptime Float: type) Float {
    comptime if (!meta.isNumeric(Float) or meta.numericType(Float) != .float)
        @compileError("zsl.float.e: Float must be a float type, got \n\tFloat = " ++ @typeName(Float) ++ "\n");

    return 2.7182818284590452353602874713526624977572;
}

// Basic functions
pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.Add: at least one of X or Y must be a float type, the other must be a bool, an int or a float type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return float.Coerce(X, Y);
}

/// Performs addition between two operands of float, int or bool types, where at
/// least one operand must be of float type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then adding them.
///
/// ## Signature
/// ```zig
/// float.add(x: X, y: Y) float.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
pub fn add(x: anytype, y: anytype) float.Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = float.Add(@TypeOf(x), @TypeOf(y));

    return numeric.cast(R, x) + numeric.cast(R, y);
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.Sub: at least one of X or Y must be a float type, the other must be a bool, an int or a float type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return float.Coerce(X, Y);
}

/// Performs subtraction between two operands of float, int or bool types, where
/// at least one operand must be of float type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then subtracting them.
///
/// ## Signature
/// ```zig
/// float.sub(x: X, y: Y) float.Sub(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
pub fn sub(x: anytype, y: anytype) float.Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = float.Sub(@TypeOf(x), @TypeOf(y));

    return numeric.cast(R, x) - numeric.cast(R, y);
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.Mul: at least one of X or Y must be a float type, the other must be a bool, an int or a float type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return float.Coerce(X, Y);
}

/// Performs multiplication between two operands of float, int or bool types,
/// where at least one operand must be of float type. The result type is
/// determined by coercing the operand types, and the operation is performed by
/// casting both operands to the result type, then multiplying them.
///
/// ## Signature
/// ```zig
/// float.mul(x: X, y: Y) float.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
pub fn mul(x: anytype, y: anytype) float.Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = float.Mul(@TypeOf(x), @TypeOf(y));

    return numeric.cast(R, x) * numeric.cast(R, y);
}

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or !meta.numericType(Z).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float and meta.numericType(Z) != .float))
        @compileError("zsl.float.Fma: at least one of X, Y or Z must be a float type, the others must be bool, int or float types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    return float.Coerce(X, numeric.Coerce(Y, Z));
}

/// Performs fused multiplication and addition (x * y + z) between three
/// operands of float, int or bool types, where at least one operand must be of
/// float type. The result type is determined by coercing the operand types, and
/// the operation is performed by casting all three operands to the result type,
/// then performing the fused operation.
///
/// ## Signature
/// ```zig
/// float.fma(x: X, y: Y, z: Z) float.Fma(X, Y, Z)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left multiplication operand.
/// * `y` (`anytype`): The right multiplication operand.
/// * `z` (`anytype`): The addition operand.
///
/// ## Returns
/// `float.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z))`: The result of the fused
/// multiplication and addition.
pub fn fma(x: anytype, y: anytype, z: anytype) float.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const R: type = float.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z));

    return @mulAdd(R, numeric.cast(R, x), numeric.cast(R, y), numeric.cast(R, z));
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.Div: at least one of X or Y must be a float type, the other must be a bool, an int or a float type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return float.Coerce(X, Y);
}

/// Performs division between two operands of float, int or bool types, where at
/// least one operand must be of float type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then dividing them.
///
/// ## Signature
/// ```zig
/// float.div(x: X, y: Y) float.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
pub fn div(x: anytype, y: anytype) float.Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = float.Div(@TypeOf(x), @TypeOf(y));

    return numeric.cast(R, x) / numeric.cast(R, y);
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for ordering. The operation is performed by
/// casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// float.order(x: X, y: Y) std.math.Order
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `std.math.Order`: The result of the comparison.
pub fn order(x: anytype, y: anytype) std.math.Order {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.order: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    if (numeric.cast(C, x) < numeric.cast(C, y)) return .lt;
    if (numeric.cast(C, x) > numeric.cast(C, y)) return .gt;
    return .eq;
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for equality. The operation is performed by
/// casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// float.eq(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are equal, `false` otherwise.
pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.eq: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    return numeric.cast(C, x) == numeric.cast(C, y);
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for inequality. The operation is performed by
/// casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// float.ne(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are not equal, `false` otherwise.
pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.ne: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    return numeric.cast(C, x) != numeric.cast(C, y);
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for for less-than ordering. The operation is
/// performed by casting both operands to the coerced type, then comparing
/// them.
///
/// ## Signature
/// ```zig
/// float.lt(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is less than `y`, `false` otherwise.
pub fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.lt: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    return numeric.cast(C, x) < numeric.cast(C, y);
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for less-than or equal ordering. The
/// operation is performed by casting both operands to the coerced type, then
/// comparing them.
///
/// ## Signature
/// ```zig
/// float.le(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is less than or equal to `y`, `false` otherwise.
pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.le: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    return numeric.cast(C, x) <= numeric.cast(C, y);
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for greater-than ordering. The operation is
/// performed by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// float.gt(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is greater than `y`, `false` otherwise.
pub fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.gt: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    return numeric.cast(C, x) > numeric.cast(C, y);
}

/// Compares two operands of float, int or bool types, where at least one
/// operand must be of float type, for greater-than or equality ordering. The
/// operation is performed by casting both operands to the coerced type, then
/// comparing them.
///
/// ## Signature
/// ```zig
/// float.ge(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if `x` is greater than or equal to `y`, `false` otherwise.
pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.ge: at least one of x or y must be a float, the other must be a bool, an int or a float, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const C: type = float.Coerce(X, Y);

    return numeric.cast(C, x) >= numeric.cast(C, y);
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.Max: at least one of X or Y must be a float type, the other must be a bool, an int or a float type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return float.Coerce(X, Y);
}

/// Returns the maximum of two operands of float, int or bool types, where at
/// least one operand must be of float type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then comparing them.
///
/// ## Signature
/// ```zig
/// float.max(x: X, y: Y) float.Max(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Max(@TypeOf(x), @TypeOf(y))`: The maximum of the two operands.
pub fn max(x: anytype, y: anytype) float.Max(@TypeOf(x), @TypeOf(y)) {
    const R: type = float.Max(@TypeOf(x), @TypeOf(y));

    return if (numeric.cast(R, x) > numeric.cast(R, y)) numeric.cast(R, x) else numeric.cast(R, y);
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.float) or !meta.numericType(Y).le(.float) or
        (meta.numericType(X) != .float and meta.numericType(Y) != .float))
        @compileError("zsl.float.Min: at least one of X or Y must be a float type, the other must be a bool, an int or a float type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return float.Coerce(X, Y);
}

/// Returns the minimum of two operands of float, int or bool types, where at
/// least one operand must be of float type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then comparing them.
///
/// ## Signature
/// ```zig
/// float.min(x: X, y: Y) float.Min(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `float.Min(@TypeOf(x), @TypeOf(y))`: The minimum of the two operands.
pub fn min(x: anytype, y: anytype) float.Min(@TypeOf(x), @TypeOf(y)) {
    const R: type = float.Min(@TypeOf(x), @TypeOf(y));

    return if (numeric.cast(R, x) < numeric.cast(R, y)) numeric.cast(R, x) else numeric.cast(R, y);
}

/// Returns the highest representable value of the given float type `Float`.
///
/// ## Arguments
/// * `Float` (`type`): The float type to get the maximum value for.
///
/// ## Returns
/// `Float`: The maximum representable value of type `Float`.
pub fn highest(comptime Float: type) Float {
    comptime if (!meta.isNumeric(Float) or meta.numericType(Float) != .float)
        @compileError("zsl.float.highest: Float must be a float type, got \n\tFloat = " ++ @typeName(Float) ++ "\n");

    return std.math.floatMax(Float);
}

/// Returns the lowest representable value of the given float type `Float`.
///
/// ## Arguments
/// * `Float` (`type`): The Float type to get the lowest value for.
///
/// ## Returns
/// `Float`: The minimum representable value of type `Float`.
pub fn lowest(comptime Float: type) Float {
    comptime if (!meta.isNumeric(Float) or meta.numericType(Float) != .float)
        @compileError("zsl.float.lowest: Float must be a float type, got \n\tFloat = " ++ @typeName(Float) ++ "\n");

    return -std.math.floatMax(Float);
}

/// Returns the smallest positive magnitude strictly greater than zero of the
/// given float type `Float`.
///
/// ## Arguments
/// * `Float` (`type`): The float type to get the minimum value for.
///
/// ## Returns
/// `Float`: The minimum representable value of type `Float`.
pub fn smallest(comptime Float: type) Float {
    comptime if (!meta.isNumeric(Float) or meta.numericType(Float) != .float)
        @compileError("zsl.float.smallest: Float must be a float type, got \n\tFloat = " ++ @typeName(Float) ++ "\n");

    return std.math.floatMin(Float);
}

pub fn isNan(x: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or meta.numericType(X) != .float)
        @compileError("zsl.float.isNan: x must be a float, got \n\tx: " ++ @typeName(X) ++ "\n");

    return std.math.isNan(x);
}

// Basic operations
pub const abs = @import("float/abs.zig").abs;
pub const sign = @import("float/sign.zig").sign;
// pub const fmod = @import("float/fmod.zig").fmod; // to implement
// pub const remainder = @import("float/remainder.zig").remainder; // to implement
// pub const remquo = @import("float/remquo.zig").remquo; // to implement
// pub const fdim = @import("float/fdim.zig").fdim; // to implement

// Exponential functions
pub const exp = @import("float/exp.zig").exp;
pub const exp2 = @import("float/exp2.zig").exp2;
pub const expm1 = @import("float/expm1.zig").expm1;
pub const ln = @import("float/ln.zig").ln;
pub const log10 = @import("float/log10.zig").log10;
pub const log2 = @import("float/log2.zig").log2;
pub const log1p = @import("float/log1p.zig").log1p;

// Power functions
pub const Pow = @import("float/pow.zig").Pow;
pub const pow = @import("float/pow.zig").pow;
pub const sqrt = @import("float/sqrt.zig").sqrt;
pub const cbrt = @import("float/cbrt.zig").cbrt;
pub const Hypot = @import("float/hypot.zig").Hypot;
pub const hypot = @import("float/hypot.zig").hypot;

// Trigonometric functions
pub const sin = @import("float/sin.zig").sin;
pub const cos = @import("float/cos.zig").cos;
pub const tan = @import("float/tan.zig").tan;
pub const asin = @import("float/asin.zig").asin;
pub const acos = @import("float/acos.zig").acos;
pub const atan = @import("float/atan.zig").atan;
pub const Atan2 = @import("float/atan2.zig").Atan2;
pub const atan2 = @import("float/atan2.zig").atan2;
pub const Sincos = @import("float/sincos.zig").Sincos;
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
