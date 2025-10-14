const types = @import("types.zig");
const scast = types.scast;
const Coerce = types.Coerce;
const Cmp = types.Cmp;

pub inline fn add(
    x: anytype,
    y: anytype,
    comptime mode: Mode,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (!(types.numericType(X) == .bool and types.numericType(X) == .int) and
        !(types.numericType(Y) == .bool and types.numericType(Y) == .int) and
        !(types.numericType(X) == .int and types.numericType(Y) == .int))
        @compileError("int.add requires at least one of x or y to be an int, the other must be a bool or an int, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime mode) {
        .default => return scast(C, x) + scast(C, y),
        .wrap => return scast(C, x) +% scast(C, y),
        .saturate => return scast(C, x) +| scast(C, y),
    }
}

pub inline fn sub(
    x: anytype,
    y: anytype,
    comptime mode: Mode,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (!(types.numericType(X) == .bool and types.numericType(X) == .int) and
        !(types.numericType(Y) == .bool and types.numericType(Y) == .int) and
        !(types.numericType(X) == .int and types.numericType(Y) == .int))
        @compileError("int.sub requires at least one of x or y to be an int, the other must be a bool or an int, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime mode) {
        .default => return scast(C, x) - scast(C, y),
        .wrap => return scast(C, x) -% scast(C, y),
        .saturate => return scast(C, x) -| scast(C, y),
    }
}

pub inline fn mul(
    x: anytype,
    y: anytype,
    comptime mode: Mode,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (!(types.numericType(X) == .bool and types.numericType(X) == .int) and
        !(types.numericType(Y) == .bool and types.numericType(Y) == .int) and
        !(types.numericType(X) == .int and types.numericType(Y) == .int))
        @compileError("int.mul requires at least one of x or y to be an int, the other must be a bool or an int, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime mode) {
        .default => return scast(C, x) * scast(C, y),
        .wrap => return scast(C, x) *% scast(C, y),
        .saturate => return scast(C, x) *| scast(C, y),
    }
}

pub inline fn div(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (!(types.numericType(X) == .bool and types.numericType(X) == .int) and
        !(types.numericType(Y) == .bool and types.numericType(Y) == .int) and
        !(types.numericType(X) == .int and types.numericType(Y) == .int))
        @compileError("int.div requires at least one of x or y to be an int, the other must be a bool or an int, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return @divTrunc(scast(C, x), scast(C, y));
}

pub inline fn cmp(
    x: anytype,
    y: anytype,
) Cmp {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.cmp requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

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

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.eq requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return x == y;
}

pub inline fn ne(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.ne requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return x != y;
}

pub inline fn lt(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.lt requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return x < y;
}

pub inline fn le(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.le requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return x <= y;
}

pub inline fn gt(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.gt requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return x > y;
}

pub inline fn ge(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.ge requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return x >= y;
}

pub inline fn max(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.max requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return if (x > y) scast(C, x) else scast(C, y);
}

pub inline fn min(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if (types.numericType(X) != .int or types.numericType(Y) != .int)
        @compileError("int.max requires both x and y to be int types, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    return if (x < y) scast(C, x) else scast(C, y);
}

pub inline fn maxVal(comptime T: type) T {
    comptime if (types.numericType(T) != .int)
        @compileError("int.max requires T to be an int type, got " ++ @typeName(T));

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return (1 << (bits - scast(@TypeOf(bits), info.int.signedness == .signed))) - 1;
}

pub inline fn minVal(comptime T: type) T {
    comptime if (types.numericType(T) != .int)
        @compileError("int.min requires T to be an int type, got " ++ @typeName(T));

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return if (info.int.signedness == .signed) -(1 << (bits - 1)) else 0;
}

pub const abs = @import("int/abs.zig").abs;

pub const Mode = enum {
    default,
    wrap,
    saturate,
};
