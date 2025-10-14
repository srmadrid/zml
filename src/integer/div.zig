const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn div(allocator: std.mem.Allocator, x: anytype, y: anytype) !Integer {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!(types.numericType(X) == .integer and types.numericType(Y) == .int) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .float) and
        !(types.numericType(X) == .integer and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .int and types.numericType(Y) == .integer) and
        !(types.numericType(X) == .float and types.numericType(Y) == .integer))
        @compileError("integer.div requires at least one of x or y to be an integer, the other must be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                var result: Integer = try .init(allocator, 0);
                errdefer result.deinit(allocator);

                try integer.div_(allocator, &result, x, y);

                return result;
            },
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return div(allocator, x, temp);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return div(allocator, temp, y);
            },
            else => unreachable,
        },
        else => unreachable,
    }
}
