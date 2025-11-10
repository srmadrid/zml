const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");
const expression = @import("../expression.zig");

/// The return type of the `eq` routine for inputs of types `X` and `Y`.
pub fn Eq(X: type, Y: type) type {
    return switch (comptime types.domainType(types.Coerce(X, Y))) {
        .expression => bool,
        .array => types.EnsureArray(types.Coerce(X, Y), bool),
        .matrix => bool,
        .vector => bool,
        .numeric => bool,
    };
}

///
pub inline fn eq(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Eq(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.eq not implemented yet for expression types"),
        .array => switch (comptime types.domainType(Y)) {
            .expression => @compileError("zml.eq not implemented yet for expression types"),
            .array, .numeric => { // array == array, array == numeric
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return array.eq(
                    ctx.array_allocator,
                    x,
                    y,
                );
            },
            else => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .expression => @compileError("zml.eq not implemented yet for expression types"),
            .array => { // numeric == array
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return array.eq(
                    ctx.array_allocator,
                    x,
                    y,
                );
            },
            .numeric => { // numeric == numeric
                switch (comptime types.numericType(C)) {
                    .bool => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return x == y;
                    },
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.eq(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.eq(x, y);
                    },
                    .cfloat => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return cfloat.eq(x, y);
                    },
                    else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
