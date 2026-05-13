//! Namespace for int operations.

const int = @This();

const options = @import("options");

const meta = @import("meta.zig");
const Cmp = meta.Cmp;
const numeric = @import("numeric.zig");

pub const Coerce = @import("int/coerce.zig").Coerce;

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.add: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return int.Coerce(X, Y);
}

/// Performs addition between two operands of int or bool types, where at least
/// one operand must be of int type. The result type is determined by coercing
/// the operand types, and the operation is performed by casting both operands
/// to the result type, then adding them.
///
/// Depending on the global int operation mode set in `options.int_mode`,
/// the addition will behave differently in case of over/underflow:
/// * `.default`: Standard addition, which will panic on over/underflow.
/// * `.wrap`: Addition with wrapping behavior on over/underflow.
/// * `.saturate`: Addition with saturation behavior on over/underflow.
///
/// ## Signature
/// ```zig
/// int.add(x: X, y: Y) int.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `int.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
pub fn add(x: anytype, y: anytype) int.Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Add(@TypeOf(x), @TypeOf(y));

    switch (comptime options.int_mode) {
        .default => return numeric.cast(R, x) + numeric.cast(R, y),
        .wrap => return numeric.cast(R, x) +% numeric.cast(R, y),
        .saturate => return numeric.cast(R, x) +| numeric.cast(R, y),
    }
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.sub: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return int.Coerce(X, Y);
}

/// Performs subtraction between two operands of int or bool types, where at
/// least one operand must be of int type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then subtracting them.
///
/// Depending on the global int operation mode set in `options.int_mode`,
/// the subtraction will behave differently in case of over/underflow:
/// * `.default`: Standard subtraction, which will panic on over/underflow.
/// * `.wrap`: Subtraction with wrapping behavior on over/underflow.
/// * `.saturate`: Subtraction with saturation behavior on over/underflow.
///
/// ## Signature
/// ```zig
/// int.sub(x: X, y: Y) int.Sub(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `int.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
pub fn sub(x: anytype, y: anytype) int.Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Sub(@TypeOf(x), @TypeOf(y));

    switch (comptime options.int_mode) {
        .default => return numeric.cast(R, x) - numeric.cast(R, y),
        .wrap => return numeric.cast(R, x) -% numeric.cast(R, y),
        .saturate => return numeric.cast(R, x) -| numeric.cast(R, y),
    }
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.mul: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return int.Coerce(X, Y);
}

/// Performs multiplication between two operands of int or bool types, where at
/// least one operand must be of int type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then multiplying them.
///
/// Depending on the global int operation mode set in `options.int_mode`,
/// the multiplication will behave differently in case of over/underflow:
/// * `.default`: Standard multiplication, which will panic on over/underflow.
/// * `.wrap`: Multiplication with wrapping behavior on over/underflow.
/// * `.saturate`: Multiplication with saturation behavior on over/underflow.
///
/// ## Signature
/// ```zig
/// int.mul(x: X, y: Y) int.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `int.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
pub fn mul(x: anytype, y: anytype) int.Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Mul(@TypeOf(x), @TypeOf(y));

    switch (comptime options.int_mode) {
        .default => return numeric.cast(R, x) * numeric.cast(R, y),
        .wrap => return numeric.cast(R, x) *% numeric.cast(R, y),
        .saturate => return numeric.cast(R, x) *| numeric.cast(R, y),
    }
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.div: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return int.Coerce(X, Y);
}

/// Performs truncating division (rounding towards zero) between two operands of
/// int or bool types, where at least one operand must be of int type. The
/// result type is determined by coercing the operand types, and the operation
/// is performed by casting both operands to the result type, then dividing
/// them.
///
/// ## Signature
/// ```zig
/// int.div(x: X, y: Y) int.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `int.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
pub fn div(x: anytype, y: anytype) int.Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Div(@TypeOf(x), @TypeOf(y));

    return @divTrunc(numeric.cast(R, x), numeric.cast(R, y));
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for ordering. The operation is performed by casting both
/// operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.cmp(x: X, y: Y) Cmp
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `Cmp`: The result of the comparison.
pub fn cmp(x: anytype, y: anytype) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.cmp: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    if (numeric.cast(C, x) < numeric.cast(C, y)) return .lt;
    if (numeric.cast(C, x) > numeric.cast(C, y)) return .gt;
    return .eq;
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for equality. The operation is performed by casting both
/// operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.eq(x: X, y: Y) bool
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
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.eq: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    return numeric.cast(C, x) == numeric.cast(C, y);
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for inequality. The operation is performed by casting both
/// operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.ne(x: X, y: Y) bool
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
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.ne: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    return numeric.cast(C, x) != numeric.cast(C, y);
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for less-than ordering. The operation is performed by
/// casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.lt(x: X, y: Y) bool
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
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.lt: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    return numeric.cast(C, x) < numeric.cast(C, y);
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for less-than or equal ordering. The operation is performed
/// by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.le(x: X, y: Y) bool
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
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.le: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    return numeric.cast(C, x) <= numeric.cast(C, y);
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for greater-than ordering. The operation is performed by
/// casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.gt(x: X, y: Y) bool
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
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.gt: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    return numeric.cast(C, x) > numeric.cast(C, y);
}

/// Compares two operands of int or bool types, where at least one operand must
/// be of int type, for greater-than or equality ordering. The operation is
/// performed by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.ge(x: X, y: Y) bool
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
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.ge: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = int.Coerce(X, Y);

    return numeric.cast(C, x) >= numeric.cast(C, y);
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.max: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return int.Coerce(X, Y);
}

/// Returns the maximum of two operands of int or bool types, where at least
/// one operand must be of int type. The result type is determined by coercing
/// the operand types, and the operation is performed by casting both operands
/// to the result type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.max(x: X, y: Y) int.Max(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `int.Max(@TypeOf(x), @TypeOf(y))`: The maximum of the two operands.
pub fn max(x: anytype, y: anytype) int.Max(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Max(@TypeOf(x), @TypeOf(y));

    return if (numeric.cast(R, x) > numeric.cast(R, y)) numeric.cast(R, x) else numeric.cast(R, y);
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.int) or !meta.numericType(Y).le(.int) or
        (meta.numericType(X) != .int and meta.numericType(Y) != .int))
        @compileError("zsl.int.min: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return int.Coerce(X, Y);
}

/// Returns the minimum of two operands of int or bool types, where at least
/// one operand must be of int type. The result type is determined by coercing
/// the operand types, and the operation is performed by casting both operands
/// to the result type, then comparing them.
///
/// ## Signature
/// ```zig
/// int.min(x: X, y: Y) int.Min(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `int.Min(@TypeOf(x), @TypeOf(y))`: The minimum of the two operands.
pub fn min(x: anytype, y: anytype) int.Min(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Min(@TypeOf(x), @TypeOf(y));

    return if (numeric.cast(R, x) < numeric.cast(R, y)) numeric.cast(R, x) else numeric.cast(R, y);
}

/// Returns the maximum representable value of the given int type `Int`.
///
/// ## Arguments
/// * `Int` (`type`): The int type to get the maximum value for.
///
/// ## Returns
/// `Int`: The maximum representable value of type `Int`.
pub fn maxVal(comptime Int: type) Int {
    comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
        @compileError("zsl.int.maxVal: Int must be an int type, got \n\tInt = " ++ @typeName(Int) ++ "\n");

    const info = @typeInfo(Int);
    const bits = info.int.bits;
    return (1 << (bits - numeric.cast(comptime_int, info.int.signedness == .signed))) - 1;
}

/// Returns the minimum representable value of the given int type `Int`.
///
/// ## Arguments
/// * `Int` (`type`): The int type to get the minimum value for.
///
/// ## Returns
/// `Int`: The minimum representable value of type `Int`.
pub fn minVal(comptime Int: type) Int {
    comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
        @compileError("zsl.int.minVal: Int must be an int type, got \n\tInt = " ++ @typeName(Int) ++ "\n");

    const info = @typeInfo(Int);
    const bits = info.int.bits;
    return if (info.int.signedness == .signed) -(1 << (bits - 1)) else 0;
}

pub const abs = @import("int/abs.zig").abs;
pub const sign = @import("int/sign.zig").sign;
pub const Pow = @import("int/pow.zig").Pow;
pub const pow = @import("int/pow.zig").pow;
