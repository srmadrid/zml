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
pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.gt_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(O)) {
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => switch (comptime types.domainType(Y)) {
                .array, .numeric => { // array = array > array, array = numeric > array, array = array > numeric, array = numeric > numeric
                    comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    } else {
                        types.validateContext(@TypeOf(ctx), .{});
                    };

                    return array.gt_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => switch (comptime types.domainType(Y)) {
                .numeric => { // numeric = numeric > numeric
                    switch (comptime types.numericType(C)) {
                        .bool => @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .int => {
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
                                int.gt(x, y),
                                ctx,
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
                                float.gt(x, y),
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
                                cfloat.gt(x, y),
                                ctx,
                            );
                        },
                        .integer => @compileError("zml.gt_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .rational => @compileError("zml.gt_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .real => @compileError("zml.gt_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => @compileError("zml.gt_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .expression => @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                    }
                },
                else => @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.gt_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
    }
}
