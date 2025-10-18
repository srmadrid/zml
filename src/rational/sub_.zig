const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const Rational = rational.Rational;

pub fn sub_(allocator: std.mem.Allocator, o: *Rational, x: anytype, y: anytype) !void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.numericType(X) != .rational and types.numericType(X) != .integer and types.numericType(X) != .int and types.numericType(X) != .float and
        types.numericType(Y) != .rational and types.numericType(Y) != .integer and types.numericType(Y) != .int and types.numericType(Y) != .float)
        @compileError("rational.sub_ requires x and y to be an int, float, integer or rational, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (!o.flags.writable)
        return rational.Error.NotWritable;

    switch (comptime types.numericType(X)) {
        .rational => switch (comptime types.numericType(Y)) {
            .rational => {
                return rational.add_(allocator, o, x, rational.neg(null, y) catch unreachable);
            },
            .integer => {
                return rational.add_(allocator, o, x, (integer.neg(null, y) catch unreachable).asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, x, rational.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .rational => {
                return rational.add_(allocator, o, x.asRational(), rational.neg(null, y) catch unreachable);
            },
            .integer => {
                return rational.add_(allocator, o, x.asRational(), (integer.neg(null, y) catch unreachable).asRational());
            },
            .float, .int => {
                var temp: Rational = try .initSet(allocator, y, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, x.asRational(), rational.neg(null, temp) catch unreachable);
            },
            else => unreachable,
        },
        .float, .int => switch (comptime types.numericType(Y)) {
            .rational => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, temp, rational.neg(null, y) catch unreachable);
            },
            .integer => {
                var temp: Rational = try .initSet(allocator, x, 1);
                defer temp.deinit(allocator);
                return rational.add_(allocator, o, temp, (integer.neg(null, y) catch unreachable).asRational());
            },
            .float, .int => {
                var tx: Rational = try .initSet(allocator, x, 1);
                defer tx.deinit(allocator);
                var ty: Rational = try .initSet(allocator, y, 1);
                defer ty.deinit(allocator);
                return rational.add_(allocator, o, tx, rational.neg(null, ty) catch unreachable);
            },
            else => unreachable,
        },

        else => unreachable,
    }
}
