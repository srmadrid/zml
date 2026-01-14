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
pub inline fn atan2_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.atan2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domain(O)) {
        .array => switch (comptime types.domain(X)) {
            .array, .numeric => switch (comptime types.domain(Y)) {
                .array, .numeric => { // array = atan2(array, array), array = atan2(numeric, array), array = atan2(array, numeric), array = atan2(numeric, numeric)
                    comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                        if (types.isArbitraryPrecision(types.Numeric(C))) {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
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
                    } else {
                        if (types.isArbitraryPrecision(types.Numeric(C))) {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        } else {
                            types.validateContext(@TypeOf(ctx), .{});
                        }
                    };

                    return array.atan2_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domain(X)) {
            .numeric => switch (comptime types.domain(Y)) {
                .numeric => { // numeric = atan2(numeric, numeric)
                    switch (comptime types.numericType(C)) {
                        .bool => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
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
                                float.atan2(x, y),
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
                                float.atan2(x, y),
                                ctx,
                            );
                        },
                        .integer => @compileError("zml.atan2_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .rational => @compileError("zml.atan2_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .real => @compileError("zml.atan2_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => @compileError("zml.atan2_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .expression => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                    }
                },
                else => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        else => @compileError("zml.atan2_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
    }
}
