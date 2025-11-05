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

/// The return type of the `sub` routine for inputs of types `X` and `Y`.
pub fn Sub(X: type, Y: type) type {
    return types.Coerce(X, Y);
}

/// Performs subtraction between two operands of compatible types.
///
/// The `sub` routine computes the subtraction `x - y`, automatically coercing
/// compatible operand types and validating the provided context. The operation
/// is performed in the coerced precision of the operands, and the resulting
/// value is returned as a new value. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric - Numeric**: scalar subtraction.
/// - **Vector - Vector**: element-wise subtraction between vectors of equal
///   length.
/// - **Matrix - Matrix**: element-wise subtraction between matrices of equal
///   shape.
/// - **Numeric - Array**, **Array - Numeric**, and **Array - Array**:
///   broadcasted element-wise subtraction.
///
/// Signature
/// ---------
/// ```zig
/// fn sub(x: X, y: Y, ctx: anytype) !Sub(X, Y)
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
/// `Sub(@TypeOf(x), @TypeOf(y))`:
/// The result of the subtraction.
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
/// `matrix.Error.DimensionMismatch`:
/// If the two matrices do not have the same shape. Can only happen if both
/// operands are matrices.
///
/// `vector.Error.DimensionMismatch`:
/// If the two vectors do not have the same length. Can only happen if both
/// operands are vectors.
pub inline fn sub(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Sub(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isMatrix(X) and !types.isMatrix(Y) and
        !types.isVector(X) and !types.isVector(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array, .numeric => { // array - array, array - numeric
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    },
                    .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex, .expression => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return array.sub(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .matrix => switch (comptime types.domainType(Y)) {
            .matrix => { // matrix - matrix
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    },
                    .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex, .expression => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return matrix.sub(
                    ctx.matrix_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"matrix_allocator"}),
                );
            },
            else => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .vector => switch (comptime types.domainType(Y)) {
            .vector => { // vector - vector
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    },
                    .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex, .expression => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return vector.sub(
                    ctx.vector_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"vector_allocator"}),
                );
            },
            else => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .array => { // numeric - array
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );
                    },
                    .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex, .expression => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return array.sub(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // numeric - numeric
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .mode = .{ .type = int.Mode, .required = false },
                            },
                        );

                        return int.sub(x, y, types.getFieldOrDefault(ctx, "mode", int.Mode, .default));
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.sub(x, y);
                    },
                    .cfloat => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return cfloat.sub(x, y);
                    },
                    .integer => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return integer.sub(ctx.allocator, x, y);
                    },
                    .rational => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return rational.sub(ctx.allocator, x, y);
                    },
                    .real => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                    .complex => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return complex.sub(ctx.allocator, x, y);
                    },
                    .expression => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
    }
}
