const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const real = @import("../real.zig");
const complex = @import("../complex.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

///
pub inline fn div_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.div_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isMatrix(O) and !types.isMatrix(X) and !types.isMatrix(Y) and
        !types.isVector(O) and !types.isVector(X) and !types.isVector(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(O, types.Coerce(X, Y));

    switch (comptime types.domainType(O)) {
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => switch (comptime types.domainType(Y)) {
                .array, .numeric => { // array = array/array, array = numeric/array, array = array/numeric, array = numeric/numeric
                    comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
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

                    return array.div_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .matrix => switch (comptime types.domainType(X)) {
            .matrix => switch (comptime types.domainType(Y)) {
                .numeric => { // matrix = matrix/numeric
                    comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
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

                    @compileError("matrix.div_ not implemented yet");
                    // return matrix.div_(
                    //     o,
                    //     x,
                    //     y,
                    //     ctx,
                    // );
                },
                else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .vector => switch (comptime types.domainType(X)) {
            .vector => switch (comptime types.domainType(Y)) {
                .numeric => { // vector = vector/numeric
                    comptime if (types.isArbitraryPrecision(types.Numeric(O))) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
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

                    @compileError("vector.div_ not implemented yet");
                    // return vector.div_(
                    //     o,
                    //     x,
                    //     y,
                    //     ctx,
                    // );
                },
                else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => switch (comptime types.domainType(Y)) {
                .numeric => { // numeric = numeric/numeric
                    switch (comptime types.numericType(C)) {
                        .bool => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .int => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            try ops.set(
                                o,
                                int.div(x, y),
                                ctx,
                            );
                        },
                        .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            try ops.set(
                                o,
                                float.div(x, y),
                                ctx,
                            );
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            try ops.set(
                                o,
                                cfloat.div(x, y),
                                ctx,
                            );
                        },
                        .integer => {
                            if (comptime types.isArbitraryPrecision(O)) {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try integer.div_(ctx.allocator, o, x, y);
                            } else {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                var result: integer.Integer = try integer.div(
                                    ctx.internal_allocator,
                                    x,
                                    y,
                                );
                                defer result.deinit(ctx.internal_allocator);

                                ops.set(o, result, .{}) catch unreachable;
                            }
                        },
                        .rational => {
                            if (comptime types.isArbitraryPrecision(O)) {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                if (comptime O == rational.Rational) {
                                    try rational.div_(ctx.allocator, o, x, y);
                                } else {
                                    var result: rational.Rational = try rational.div(
                                        ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(ctx.allocator);

                                    try ops.set(o, result, .{ .allocator = ctx.allocator });
                                }
                            } else {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                var result: rational.Rational = try rational.div(
                                    ctx.internal_allocator,
                                    x,
                                    y,
                                );
                                defer result.deinit(ctx.internal_allocator);

                                ops.set(o, result, .{}) catch unreachable;
                            }
                        },
                        .real => @compileError("zml.div_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => {
                            if (comptime types.isArbitraryPrecision(O)) {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                if (comptime O == complex.Complex(rational.Rational) or O == complex.Complex(real.Real)) {
                                    try complex.div_(ctx.allocator, o, x, y);
                                } else {
                                    var result: types.Coerce(X, Y) = try complex.div(
                                        ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(ctx.allocator);

                                    try ops.set(o, result, .{ .allocator = ctx.allocator });
                                }
                            } else {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                var result: types.Coerce(X, Y) = try complex.div(
                                    ctx.internal_allocator,
                                    x,
                                    y,
                                );
                                defer result.deinit(ctx.internal_allocator);

                                ops.set(o, result, .{}) catch unreachable;
                            }
                        },
                        .expression => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                    }
                },
                else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.div_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
    }
}
