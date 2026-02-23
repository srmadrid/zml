//! Namespace for int operations.

const int = @This();

const options = @import("options");

const types = @import("types.zig");
const Cmp = types.Cmp;

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.add: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

/// Performs addition between two operands of int or bool types, where at least
/// one operand must be of int type. The result type is determined by coercing
/// the operand types, and the operation is performed by casting both operands
/// to the result type, then adding them.
///
/// Depending on the global int operation mode set in `options.int_mode`,
/// the addition will behave differently in case of over/underflow:
/// - `.default`: Standard addition, which will panic on over/underflow.
/// - `.wrap`: Addition with wrapping behavior on over/underflow.
/// - `.saturate`: Addition with saturation behavior on over/underflow.
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
pub inline fn add(x: anytype, y: anytype) int.Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Add(@TypeOf(x), @TypeOf(y));

    switch (comptime options.int_mode) {
        .default => return types.scast(R, x) + types.scast(R, y),
        .wrap => return types.scast(R, x) +% types.scast(R, y),
        .saturate => return types.scast(R, x) +| types.scast(R, y),
    }
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.sub: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

/// Performs subtraction between two operands of int or bool types, where at
/// least one operand must be of int type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then subtracting them.
///
/// Depending on the global int operation mode set in `options.int_mode`,
/// the subtraction will behave differently in case of over/underflow:
/// - `.default`: Standard subtraction, which will panic on over/underflow.
/// - `.wrap`: Subtraction with wrapping behavior on over/underflow.
/// - `.saturate`: Subtraction with saturation behavior on over/underflow.
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
pub inline fn sub(x: anytype, y: anytype) int.Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Sub(@TypeOf(x), @TypeOf(y));

    switch (comptime options.int_mode) {
        .default => return types.scast(R, x) - types.scast(R, y),
        .wrap => return types.scast(R, x) -% types.scast(R, y),
        .saturate => return types.scast(R, x) -| types.scast(R, y),
    }
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.mul: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return types.Coerce(X, Y);
}

/// Performs multiplication between two operands of int or bool types, where at
/// least one operand must be of int type. The result type is determined by
/// coercing the operand types, and the operation is performed by casting both
/// operands to the result type, then multiplying them.
///
/// Depending on the global int operation mode set in `options.int_mode`,
/// the multiplication will behave differently in case of over/underflow:
/// - `.default`: Standard multiplication, which will panic on over/underflow.
/// - `.wrap`: Multiplication with wrapping behavior on over/underflow.
/// - `.saturate`: Multiplication with saturation behavior on over/underflow.
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
pub inline fn mul(x: anytype, y: anytype) int.Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Mul(@TypeOf(x), @TypeOf(y));

    switch (comptime options.int_mode) {
        .default => return types.scast(R, x) * types.scast(R, y),
        .wrap => return types.scast(R, x) *% types.scast(R, y),
        .saturate => return types.scast(R, x) *| types.scast(R, y),
    }
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.div: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return types.Coerce(X, Y);
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
pub inline fn div(x: anytype, y: anytype) int.Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Div(@TypeOf(x), @TypeOf(y));

    return @divTrunc(types.scast(R, x), types.scast(R, y));
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
pub inline fn cmp(x: anytype, y: anytype) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.cmp: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    if (types.scast(C, x) < types.scast(C, y)) return .lt;
    if (types.scast(C, x) > types.scast(C, y)) return .gt;
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
pub inline fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.eq: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) == types.scast(C, y);
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
pub inline fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.ne: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) != types.scast(C, y);
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
pub inline fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.lt: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) < types.scast(C, y);
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
pub inline fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.le: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) <= types.scast(C, y);
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
pub inline fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.gt: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) > types.scast(C, y);
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
pub inline fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.ge: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    return types.scast(C, x) >= types.scast(C, y);
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.max: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return types.Coerce(X, Y);
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
pub inline fn max(x: anytype, y: anytype) int.Max(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Max(@TypeOf(x), @TypeOf(y));

    return if (types.scast(R, x) > types.scast(R, y)) types.scast(R, x) else types.scast(R, y);
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.min: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return types.Coerce(X, Y);
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
pub inline fn min(x: anytype, y: anytype) int.Min(@TypeOf(x), @TypeOf(y)) {
    const R: type = int.Min(@TypeOf(x), @TypeOf(y));

    return if (types.scast(R, x) < types.scast(R, y)) types.scast(R, x) else types.scast(R, y);
}

/// Returns the maximum representable value of the given int type `T`.
///
/// ## Arguments
/// * `T` (`type`): The int type to get the maximum value for.
///
/// ## Returns
/// `T`: The maximum representable value of type `T`.
pub inline fn maxVal(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .int)
        @compileError("zml.int.maxVal: T must be an int type, got \n\tT: " ++ @typeName(T) ++ "\n");

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return (1 << (bits - types.scast(@TypeOf(bits), info.int.signedness == .signed))) - 1;
}

/// Returns the minimum representable value of the given int type `T`.
///
/// ## Arguments
/// * `T` (`type`): The int type to get the minimum value for.
///
/// ## Returns
/// `T`: The minimum representable value of type `T`.
pub inline fn minVal(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .int)
        @compileError("zml.int.minVal: T must be an int type, got \n\tT: " ++ @typeName(T) ++ "\n");

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return if (info.int.signedness == .signed) -(1 << (bits - 1)) else 0;
}

pub const abs = @import("int/abs.zig").abs;
pub const sign = @import("int/sign.zig").sign;
pub const Pow = @import("int/pow.zig").Pow;
pub const pow = @import("int/pow.zig").pow;

pub const Error = error{
    NegativeExponent,
};
