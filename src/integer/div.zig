const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn div(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("integer.div requires x and y to be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    var result: Integer = try .init(allocator, 0);
    errdefer result.deinit(allocator);

    try integer.div_(allocator, &result, x, y);

    return result;
}
