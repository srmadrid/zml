const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

/// The return type of the `lt` routine for inputs of types `X` and `Y`.
pub fn Lt(X: type, Y: type) type {
    return switch (comptime types.domainType(types.Coerce(X, Y))) {
        .array => types.EnsureArray(types.Coerce(X, Y), bool),
        .matrix => bool,
        .vector => bool,
        .numeric => bool,
    };
}

///
pub inline fn lt(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Lt(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array, .numeric => { // array < array, array < numeric
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return array.lt(
                    ctx.array_allocator,
                    x,
                    y,
                );
            },
            else => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .array => { // numeric < array
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return array.lt(
                    ctx.array_allocator,
                    x,
                    y,
                );
            },
            .numeric => { // numeric < numeric
                switch (comptime types.numericType(C)) {
                    .bool => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return x < y;
                    },
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.lt(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.lt(x, y);
                    },
                    .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
