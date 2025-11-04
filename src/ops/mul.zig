const std = @import("std");

const types = @import("../types.zig");
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

/// The return type of the `mul` routine for inputs of types `X` and `Y`.
pub fn Mul(X: type, Y: type) type {
    return types.MulCoerce(X, Y);
}

/// Performs multiplication between two operands of compatible types.
///
/// The `mul` routine computes the product `x * y`, automatically coercing
/// compatible operand types and validating the provided context. The operation
/// is performed in the coerced precision of the operands, and the resulting
/// value is returned as a new value. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric * Numeric**: scalar multiplication.
/// - **Numeric * Vector** and **Vector * Numeric**: element-wise multiplication
///   of a vector by a scalar.
/// - **Numeric * Matrix** and **Matrix * Numeric**: element-wise multiplication
///   of a matrix by a scalar.
/// - **Vector * Vector**: vector dot product.
/// - **Vector * Matrix** and **Matrix * Vector**: matrix-vector multiplication,
///   with the vector treated as a row or column vector, respectively.
/// - **Matrix * Matrix**: matrix-matrix multiplication.
/// - **Numeric * Array**, **Array * Numeric**, and **Array * Array**:
///   broadcasted element-wise multiplication.
///
/// Signature
/// ---------
/// ```zig
/// fn mul(x: X, y: Y, ctx: anytype) !Mul(X, Y)
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
/// `Mul(@TypeOf(x), @TypeOf(y))`:
/// The result of the multiplication.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the coerced type is of
/// arbitrary precision or a structured data type.
///
/// `array.Error.NotBroadcastable`:
/// If the two arrays cannot be broadcasted to a common shape. Can only happen
/// if at least one of the operands is an array.
///
/// `linalg.Error.DimensionMismatch`:
/// If the dimensions of the operands are incompatible for matrix,
/// matrix-vector, vector-matrix, or vector-vector multiplication. Can only
/// happen in the above cases.
pub inline fn mul(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Mul(@TypeOf(x), @TypeOf(y)) {
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
                    .rational => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return rational.mul(ctx.allocator, x, y);
                    },
                    .real => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                    .complex => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return complex.mul(ctx.allocator, x, y);
                    },
                    .expression => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
        },
    }
}
