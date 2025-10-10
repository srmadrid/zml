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
pub inline fn min_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.min_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isMatrix(O) and !types.isMatrix(X) and !types.isMatrix(Y) and
        !types.isVector(O) and !types.isVector(X) and !types.isVector(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(O)) {
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => switch (comptime types.domainType(Y)) {
                .array, .numeric => { // array = min(array, array), array = min(numeric, array), array = min(array, numeric), array = min(numeric, numeric)
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

                    return array.min_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => switch (comptime types.domainType(Y)) {
                .numeric => { // numeric = min(numeric, numeric)
                    switch (comptime types.numericType(C)) {
                        .bool => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
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
                                int.min(x, y),
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
                                float.min(x, y),
                                ctx,
                            );
                        },
                        .cfloat => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .integer => { // Sketch. Not implemented yet.
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
                                integer.min(null, x, y),
                                ctx,
                            );
                        },
                        .rational => @compileError("zml.min_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .real => @compileError("zml.min_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .expression => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                    }
                },
                else => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        else => @compileError("zml.min_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
    }
}
