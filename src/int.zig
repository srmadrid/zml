const types = @import("types.zig");
const cast = types.cast;
const Coerce = types.Coerce;

pub fn add(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (types.numericType(@TypeOf(left)) != .int)
        @compileError("left must be an int");

    comptime if (types.numericType(@TypeOf(right)) != .int)
        @compileError("right must be an int");

    const R: type = Coerce(@TypeOf(left), @TypeOf(right));

    return cast(R, left, .{}) + cast(R, right, .{});
}

pub fn sub(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (types.numericType(@TypeOf(left)) != .int)
        @compileError("left must be an int");

    comptime if (types.numericType(@TypeOf(right)) != .int)
        @compileError("right must be an int");

    const R: type = Coerce(@TypeOf(left), @TypeOf(right));

    return cast(R, left, .{}) - cast(R, right, .{});
}

pub fn mul(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (types.numericType(@TypeOf(left)) != .int)
        @compileError("left must be an int");

    comptime if (types.numericType(@TypeOf(right)) != .int)
        @compileError("right must be an int");

    const R: type = Coerce(@TypeOf(left), @TypeOf(right));

    return cast(R, left, .{}) * cast(R, right, .{});
}

pub fn div(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (types.numericType(@TypeOf(left)) != .int)
        @compileError("left must be an int");

    comptime if (types.numericType(@TypeOf(right)) != .int)
        @compileError("right must be an int");

    const R: type = Coerce(@TypeOf(left), @TypeOf(right));

    return cast(R, left, .{}) / cast(R, right, .{});
}

pub const abs = @import("int/abs.zig").abs;
