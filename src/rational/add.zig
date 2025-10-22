const std = @import("std");

const types = @import("../types.zig");
const float = @import("../float.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

pub fn add(allocator: std.mem.Allocator, x: anytype, y: anytype) !Rational {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.add requires x and y to be an int, float, integer or rational, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    var result: Rational = try .init(allocator, 0, 0);
    errdefer result.deinit(allocator);

    try rational.add_(allocator, &result, x, y);

    return result;
}
