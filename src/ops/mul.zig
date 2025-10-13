const std = @import("std");

const types = @import("../types.zig");
const MulCoerce = types.MulCoerce;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

///
pub inline fn mul(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !MulCoerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isMatrix(X) and !types.isMatrix(Y) and
        !types.isVector(X) and !types.isVector(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.MulCoerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array, .numeric => { // array * array, array * numeric
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    if (types.numericType(types.Numeric(C)) == .int) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                };

                return array.mul(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .matrix => switch (comptime types.domainType(Y)) {
            .matrix => { // matrix * matrix
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return matrix.mul(
                    ctx.matrix_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"matrix_allocator"}),
                );
            },
            .vector => { // matrix * vector
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return matrix.mul(
                    ctx.vector_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"vector_allocator"}),
                );
            },
            .numeric => { // matrix * numeric
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    if (types.numericType(types.Numeric(C)) == .int) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                };

                return matrix.mul(
                    ctx.matrix_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"matrix_allocator"}),
                );
            },
            else => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .vector => switch (comptime types.domainType(Y)) {
            .matrix => { // vector * matrix
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return matrix.mul(
                    ctx.vector_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"vector_allocator"}),
                );
            },
            .vector => { // vector * vector
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(@TypeOf(ctx), .{});
                };

                return vector.mul(
                    types.useless_allocator,
                    x,
                    y,
                    ctx,
                );
            },
            else => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .array => { // numeric * array
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    if (types.numericType(types.Numeric(C)) == .int) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                };

                return array.mul(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .matrix => { // numeric * matrix
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    if (types.numericType(types.Numeric(C)) == .int) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                };

                return matrix.mul(
                    ctx.matrix_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"matrix_allocator"}),
                );
            },
            .vector => { // numeric * vector
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    if (types.numericType(types.Numeric(C)) == .int) {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    } else {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    }
                };

                return vector.mul(
                    ctx.vector_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"vector_allocator"}),
                );
            },
            .numeric => { // numeric * numeric
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );

                        return int.mul(x, y, types.getFieldOrDefault(ctx, "mode", int.Mode, .default));
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.mul(x, y);
                    },
                    .cfloat => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return cfloat.mul(x, y);
                    },
                    .integer => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return integer.mul(ctx.allocator, x, y);
                    },
                    else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
        },
    }
}
