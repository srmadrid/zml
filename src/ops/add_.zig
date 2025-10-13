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

/// Performs in-place element-wise addition between two operands of compatible
/// types.
///
/// The `add_` routine computes the sum `x + y` and stores the result directly
/// into `o`, automatically validating the provided context. The operation is
/// performed in the coerced precision of the operands and the output, and the
/// result is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = Numeric + Numeric**: scalar addition.
/// - **Vector = Vector + Vector**: element-wise addition between vectors of
///   equal length.
/// - **Matrix = Matrix + Matrix**: element-wise addition between matrices of
///   equal shape.
/// - **Array = Numeric + Numeric**, **Array = Numeric + Array**, **Array =
///   Array + Numeric**, and **Array = Array + Array**: broadcasted element-wise
///   addition.
///
/// Signature
/// ---------
/// ```zig
/// fn add_(o: *O, x: X, y: Y, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// `ctx` (`anytype`):
/// A context struct providing necessary resources and configuration for the
/// operation. The required fields depend on the output and operand types. If
/// the context is missing required fields or contains unnecessary or wrongly
/// typed fields, the compiler will emit a detailed error message describing the
/// expected structure.
///
/// Returns
/// -------
/// `void`:
/// The result is written in place to `o`.
///
/// Errors
/// ------
/// ``
pub inline fn add_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isMatrix(O) and !types.isMatrix(X) and !types.isMatrix(Y) and
        !types.isVector(O) and !types.isVector(X) and !types.isVector(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(O, types.Coerce(X, Y));

    switch (comptime types.domainType(O)) {
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => switch (comptime types.domainType(Y)) {
                .array, .numeric => { // array = array + array, array = numeric + array, array = array + numeric, array = numeric + numeric
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

                    return array.add_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .matrix => switch (comptime types.domainType(X)) {
            .matrix => switch (comptime types.domainType(Y)) {
                .matrix => { // matrix = matrix + matrix
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

                    @compileError("matrix.add_ not implemented yet");
                    // return matrix.add_(
                    //     o,
                    //     x,
                    //     y,
                    //     ctx,
                    // );
                },
                else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .vector => switch (comptime types.domainType(X)) {
            .vector => switch (comptime types.domainType(Y)) {
                .vector => { // vector = vector + vector
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

                    @compileError("vector.add_ not implemented yet");
                    // return vector.add_(
                    //     o,
                    //     x,
                    //     y,
                    //     ctx,
                    // );
                },
                else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => switch (comptime types.domainType(Y)) {
                .numeric => { // numeric = numeric + numeric
                    switch (comptime types.numericType(C)) {
                        .bool => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .mode = .{ .type = int.Mode, .required = false },
                                },
                            );

                            try ops.set(
                                o,
                                int.add(
                                    x,
                                    y,
                                    types.getFieldOrDefault(ctx, "mode", int.Mode, .default),
                                ),
                                types.stripStruct(ctx, &.{"mode"}),
                            );
                        },
                        .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            try ops.set(
                                o,
                                float.add(x, y),
                                ctx,
                            );
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            try ops.set(
                                o,
                                cfloat.add(x, y),
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

                                try integer.add_(ctx.allocator, o, x, y);
                            } else {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                var result: integer.Integer = try integer.add(
                                    ctx.internal_allocator,
                                    x,
                                    y,
                                );
                                defer result.deinit(ctx.internal_allocator);

                                ops.set(o, result, .{}) catch unreachable;
                            }
                        },
                        .rational => @compileError("zml.add_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .real => @compileError("zml.add_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => @compileError("zml.add_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .expression => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                    }
                },
                else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.add_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
    }
}
