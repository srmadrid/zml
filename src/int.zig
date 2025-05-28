const types = @import("types.zig");
const scast = types.scast;
const Coerce = types.Coerce;

pub inline fn add(
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: Mode = .default,
    },
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.add requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (options.mode) {
        .default => return scast(C, left) + scast(C, right),
        .wrap => return scast(C, left) +% scast(C, right),
        .saturate => return scast(C, left) +| scast(C, right),
    }
}

pub inline fn add_(
    out: anytype,
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: Mode = .default,
    },
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("int.add_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.add_ requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (options.mode) {
        .default => out.* = scast(O, scast(C, left) + scast(C, right)),
        .wrap => out.* = scast(O, scast(C, left) +% scast(C, right)),
        .saturate => out.* = scast(O, scast(C, left) +| scast(C, right)),
    }
}

pub inline fn sub(
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: Mode = .default,
    },
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.sub requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (options.mode) {
        .default => return scast(C, left) - scast(C, right),
        .wrap => return scast(C, left) -% scast(C, right),
        .saturate => return scast(C, left) -| scast(C, right),
    }
}

pub inline fn sub_(
    out: anytype,
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: Mode = .default,
    },
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("int.sub_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.sub_ requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (options.mode) {
        .default => out.* = scast(O, scast(C, left) - scast(C, right)),
        .wrap => out.* = scast(O, scast(C, left) -% scast(C, right)),
        .saturate => out.* = scast(O, scast(C, left) -| scast(C, right)),
    }
}

pub inline fn mul(
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: Mode = .default,
    },
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.mul requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (options.mode) {
        .default => return scast(C, left) * scast(C, right),
        .wrap => return scast(C, left) *% scast(C, right),
        .saturate => return scast(C, left) *| scast(C, right),
    }
}

pub inline fn mul_(
    out: anytype,
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: Mode = .default,
    },
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("int.mul_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.mul_ requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (options.mode) {
        .default => out.* = scast(O, scast(C, left) * scast(C, right)),
        .wrap => out.* = scast(O, scast(C, left) *% scast(C, right)),
        .saturate => out.* = scast(O, scast(C, left) *| scast(C, right)),
    }
}

pub inline fn div(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.div requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
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
        @compileError("int.div_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int) or
        (types.numericType(R) != .bool and types.numericType(R) != .int) or
        (types.numericType(L) != .int and types.numericType(R) != .int))
        @compileError("int.div_ requires at least one of left or right to be an int, the other must be a bool or an int, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    out.* = scast(O, scast(C, left) / scast(C, right));
}

pub inline fn max(comptime T: type) T {
    comptime if (types.numericType(T) != .int)
        @compileError("int.max requires T to be an int type, got " ++ @typeName(T));

    const info = @typeInfo(T);
    const bits = info.int.bits;
    return (1 << (bits - scast(@TypeOf(bits), info.int.signedness == .signed))) - 1;
}

pub inline fn min(comptime T: type) T {
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
