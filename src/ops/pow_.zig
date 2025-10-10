const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

///
pub inline fn pow_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.pow_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(O)) {
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => switch (comptime types.domainType(Y)) {
                .array, .numeric => { // array = array^array, array = numeric^array, array = array^numeric, array = numeric^numeric
                    comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                        if (types.isArbitraryPrecision(types.Numeric(C))) {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            if (types.numericType(types.Numeric(C)) == .int) {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .mode = .{ .type = int.Mode, .required = false },
                                    },
                                );
                            } else {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );
                            }
                        }
                    } else {
                        if (types.isArbitraryPrecision(types.Numeric(C))) {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            if (types.numericType(types.Numeric(C)) == .int) {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .mode = .{ .type = int.Mode, .required = false },
                                    },
                                );
                            } else {
                                types.validateContext(@TypeOf(ctx), .{});
                            }
                        }
                    };

                    return array.pow_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => switch (comptime types.domainType(Y)) {
                .numeric => { // numeric = numeric^numeric
                    switch (comptime types.numericType(C)) {
                        .bool => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .int => {
                            comptime if (types.isArbitraryPrecision(O)) {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .mode = .{ .type = int.Mode, .required = false },
                                    },
                                );
                            } else {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .mode = .{ .type = int.Mode, .required = false },
                                    },
                                );
                            };

                            try ops.set(
                                o,
                                int.pow(
                                    x,
                                    y,
                                    types.getFieldOrDefault(ctx, "mode", int.Mode, .default),
                                ),
                                types.stripStruct(ctx, &.{"mode"}),
                            );
                        },
                        .float => {
                            comptime if (types.isArbitraryPrecision(O)) {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );
                            } else {
                                types.validateContext(@TypeOf(ctx), .{});
                            };

                            try ops.set(
                                o,
                                float.pow(x, y),
                                ctx,
                            );
                        },
                        .cfloat => {
                            comptime if (types.isArbitraryPrecision(O)) {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );
                            } else {
                                types.validateContext(@TypeOf(ctx), .{});
                            };

                            try ops.set(
                                o,
                                cfloat.pow(x, y),
                                ctx,
                            );
                        },
                        .integer => @compileError("zml.pow_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .rational => @compileError("zml.pow_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .real => @compileError("zml.pow_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => @compileError("zml.pow_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .expression => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                    }
                },
                else => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        else => @compileError("zml.pow_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
    }
}
