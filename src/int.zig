//! Namespace for int operations.

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

pub inline fn add(x: anytype, y: anytype) Add(@TypeOf(x), @TypeOf(y)) {
    const R: type = Add(@TypeOf(x), @TypeOf(y));

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

pub inline fn sub(x: anytype, y: anytype) Sub(@TypeOf(x), @TypeOf(y)) {
    const R: type = Sub(@TypeOf(x), @TypeOf(y));

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

pub inline fn mul(x: anytype, y: anytype) Mul(@TypeOf(x), @TypeOf(y)) {
    const R: type = Mul(@TypeOf(x), @TypeOf(y));

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

pub inline fn div(x: anytype, y: anytype) Div(@TypeOf(x), @TypeOf(y)) {
    const R: type = Div(@TypeOf(x), @TypeOf(y));

    return @divTrunc(types.scast(R, x), types.scast(R, y));
}

pub inline fn cmp(
    x: anytype,
    y: anytype,
) Cmp {
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

pub inline fn eq(
    x: anytype,
    y: anytype,
) bool {
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

pub inline fn ne(
    x: anytype,
    y: anytype,
) bool {
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

pub inline fn lt(
    x: anytype,
    y: anytype,
) bool {
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

pub inline fn le(
    x: anytype,
    y: anytype,
) bool {
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

pub inline fn gt(
    x: anytype,
    y: anytype,
) bool {
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

pub inline fn ge(
    x: anytype,
    y: anytype,
) bool {
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

pub fn Max(
    comptime X: type,
    comptime Y: type,
) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.max: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return types.Coerce(X, Y);
}

pub inline fn max(
    x: anytype,
    y: anytype,
) Max(@TypeOf(x), @TypeOf(y)) {
    const R: type = Max(@TypeOf(x), @TypeOf(y));

    return if (types.scast(R, x) > types.scast(R, y)) types.scast(R, x) else types.scast(R, y);
}

pub fn Min(
    comptime X: type,
    comptime Y: type,
) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.int) or !types.numericType(Y).le(.int) or
        (types.numericType(X) != .int and types.numericType(Y) != .int))
        @compileError("zml.int.min: at least one of x or y must be an int, the other must be a bool or an int, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y));

    return types.Coerce(X, Y);
}

pub inline fn min(
    x: anytype,
    y: anytype,
) Min(@TypeOf(x), @TypeOf(y)) {
    const R: type = Min(@TypeOf(x), @TypeOf(y));

    return if (types.scast(R, x) < types.scast(R, y)) types.scast(R, x) else types.scast(R, y);
}

pub inline fn maxVal(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .int)
        @compileError("zml.int.maxVal: T must be an int type, got \n\tT: " ++ @typeName(T) ++ "\n");

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return (1 << (bits - types.scast(@TypeOf(bits), info.int.signedness == .signed))) - 1;
}

pub inline fn minVal(comptime T: type) T {
    comptime if (!types.isNumeric(T) or types.numericType(T) != .int)
        @compileError("zml.int.minVal: T must be an int type, got \n\tT: " ++ @typeName(T) ++ "\n");

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return if (info.int.signedness == .signed) -(1 << (bits - 1)) else 0;
}

pub const abs = @import("int/abs.zig").abs;
pub const Pow = @import("int/pow.zig").Pow;
pub const pow = @import("int/pow.zig").pow;

pub const Error = error{
    NegativeExponent,
};
