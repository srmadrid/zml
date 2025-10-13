const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;

pub fn sub_(allocator: std.mem.Allocator, o: *Integer, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("integer.sub_ requires x and y to be an int, float or integer, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (comptime types.numericType(X)) {
        .integer => switch (comptime types.numericType(Y)) {
            .integer => {
                return integer.add_(allocator, o, x, integer.neg(null, y) catch unreachable);
            },
            .float, .int => {
                var temp: Integer = try .initSet(allocator, y);
                defer temp.deinit(allocator);
                return integer.add_(allocator, o, x, integer.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .float, .int => {
                var tx: Integer = try .initSet(allocator, x);
                defer tx.deinit(allocator);
                var ty: Integer = try .initSet(allocator, y);
                defer ty.deinit(allocator);
                return integer.add_(allocator, o, tx, integer.neg(null, ty) catch unreachable);
            },
            .integer => {
                var temp: Integer = try .initSet(allocator, x);
                defer temp.deinit(allocator);
                return integer.add_(allocator, o, temp, integer.neg(null, y) catch unreachable);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
