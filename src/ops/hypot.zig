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

/// The return type of the `hypot` routine for inputs of types `X` and `Y`.
pub fn Hypot(X: type, Y: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => switch (comptime types.domain(Y)) {
            .expression => expression.Expression,
            .array => types.EnsureArray(Y, Hypot(types.Numeric(X), types.Numeric(Y))),
            .matrix => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .vector => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.EnsureArray(Y, Hypot(types.Numeric(X), Y)),
        },
        .matrix => switch (comptime types.domain(Y)) {
            .expression => expression.Expression,
            .array => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .matrix => @compileError("zml.Hypot not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet"),
            .vector => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .vector => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .numeric => switch (comptime types.domain(Y)) {
            .array => types.EnsureArray(Y, Hypot(X, types.Numeric(Y))),
            .matrix => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .vector => @compileError("zml.Hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.EnsureFloat(types.Coerce(X, Y)),
        },
    };
}

/// Performs the hypotenuse operation `√(x² + y²)`.
///
/// The `hypot` routine computes the hypotenuse operation `√(x² + y²)`,
/// automatically coercing compatible operand types and validating the provided
/// context. The operation is performed in the coerced precision of the
/// operands, and the resulting value is returned as a new value. It supports
/// both fixed-precision and arbitrary-precision arithmetic, as well as
/// structured data domains. The supported type combinations are:
/// - **(Numeric, Numeric)**: scalar hypotenuse.
/// - **(Matrix, Matrix)**: matrix hypotenuse (not implemented yet).
/// - **(Numeric, Array)**, **(Array, Numeric)**, and **(Array, Array)**:
///   broadcasted element-wise hypotenuse.
/// - **(Numeric, Expression)**, **(Matrix, Expression)**,
///   **(Array, Expression)**, **(Expression, Numeric)**,
///   **(Expression, Matrix)**, **(Expression, Array)**, and
///   **(Expression, Expression)**: symbolic hypotenuse.
///
/// Signature
/// ---------
/// ```zig
/// fn hypot(x: X, y: Y, ctx: anytype) !Hypot(X, Y)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// `ctx` (`anytype`):
/// A context struct providing necessary resources and configuration for the
/// operation. The required fields depend on the operand types. If the context
/// is missing required fields or contains unnecessary or wrongly typed fields,
/// the compiler will emit a detailed error message describing the expected
/// structure.
///
/// Returns
/// -------
/// `Hypot(@TypeOf(x), @TypeOf(y))`:
/// The result of the hypotenuse operation.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the coerced type is of
/// arbitrary precision or a structured data type.
///
/// `array.Error.NotBroadcastable`:
/// If the two arrays cannot be broadcasted to a common shape. Can only happen
/// if both operands are arrays.
///
/// `matrix.Error....`:
/// Tbd.
pub inline fn hypot(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Hypot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.hypot not implemented for expressions yet"),
        .array => switch (comptime types.domain(Y)) {
            .expression => @compileError("zml.hypot not implemented for expressions yet"),
            .array, .numeric => { // hypot(array, array), hypot(array, numeric)
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return array.hypot(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domain(Y)) {
            .expression => @compileError("zml.hypot not implemented for expressions yet"),
            .array => { // hypot(numeric, array)
                comptime if (types.isArbitraryPrecision(types.Numeric(C))) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                };

                return array.hypot(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // hypot(numeric, numeric)
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.hypot(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.hypot(x, y);
                    },
                    else => @compileError("zml.hypot between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.hypot not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
