const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

/// The return type of the `min` routine for inputs of types `X` and `Y`.
pub fn Min(X: type, Y: type) type {
    return types.Coerce(X, Y);
}

///
pub inline fn min(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Min(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array, .numeric => { // min(array, array), min(array, numeric)
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

                return array.min(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .array => { // min(numeric, array)
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

                return array.min(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // min(numeric, numeric)
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.min(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.min(x, y);
                    },
                    .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .integer => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                            },
                        );

                        @compileError("Not yet implemented"); // Just to show how the context might be
                        // Both are Integer, or one is an int -> can also create a view of it if the int is the min
                        //return integer.min(ctx.allocator, x, y);
                    },
                    .rational => {
                        comptime if (types.numericType(X) == .float or types.numericType(Y) == .float) {
                            // Cannot create a Rational as a view of a float, if that is the minimum. So always allocate
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
                        //return integer.min(ctx.allocator, x, y);
                    },
                    else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
