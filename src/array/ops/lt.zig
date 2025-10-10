const std = @import("std");

const types = @import("../../types.zig");
const EnsureArray = types.EnsureArray;
const Coerce = types.Coerce;
const ops = @import("../../ops.zig");

const arrops = @import("../ops.zig");

///
pub inline fn lt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y))
        @compileError("zml.array.lt requires at least one array argument");

    return arrops.apply2(
        allocator,
        x,
        y,
        ops.lt,
        .{},
    );
}
