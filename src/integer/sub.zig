const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn sub(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((types.numericType(X) != .integer and types.numericType(X) != .int) or
        (types.numericType(X) != .integer and types.numericType(X) != .float) or
        (types.numericType(X) != .int and types.numericType(X) != .integer) or
        (types.numericType(X) != .float and types.numericType(X) != .integer))
        @compileError("integer.sub requires at least one of x or y to be an integer, the other must be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return integer.add(allocator, x, ops.neg(y, .{}) catch unreachable);
}
