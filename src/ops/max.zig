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

/// The return type of the `max` routine for inputs of types `X` and `Y`.
pub fn Max(X: type, Y: type) type {
    return types.Coerce(X, Y);
}

///
pub inline fn max(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Max(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.max not implemented for expression types yet"),
        .array => switch (comptime types.domainType(Y)) {
            .expression => @compileError("zml.max not implemented for expression types yet"),
            .array, .numeric => { // max(array, array), max(array, numeric)
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    if (types.numericType(C) == .rational and (types.numericType(X) == .float or types.numericType(Y) == .float)) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return array.max(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .expression => @compileError("zml.max not implemented for expression types yet"),
            .array => { // max(numeric, array)
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    if (types.numericType(C) == .rational and (types.numericType(X) == .float or types.numericType(Y) == .float)) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return array.max(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // max(numeric, numeric)
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.max(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.max(x, y);
                    },
                    .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .integer => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                            },
                        );

                        @compileError("Not yet implemented"); // Just to show how the context might be
                        // Both are Integer, or one is an int -> can also create a view of it if the int is the max
                        //return integer.max(ctx.allocator, x, y);
                    },
                    .rational => {
                        comptime if (types.numericType(X) == .float or types.numericType(Y) == .float) {
                            // Cannot create a Rational as a view of a float, if that is the maximum. So always allocate
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            // Can create a Rational as a view of an int or Integer
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                                },
                            );
                        };

                        @compileError("Not yet implemented"); // Just to show how the context might be
                        //return integer.max(ctx.allocator, x, y);
                    },
                    else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
