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
const expression = @import("../expression.zig");

/// The return type of the `add` routine for inputs of types `X` and `Y`.
pub fn Add(X: type, Y: type) type {
    return types.Coerce(X, Y);
}

/// Performs addition between two operands of compatible types.
///
/// The `add` routine computes the sum `x + y`, automatically coercing
/// compatible operand types and validating the provided context. The operation
/// is performed in the coerced precision of the operands, and the resulting
/// value is returned as a new value. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric + Numeric**: scalar addition.
/// - **Vector + Vector**: element-wise addition between vectors of equal
///   length.
/// - **Matrix + Matrix**: element-wise addition between matrices of equal
///   shape.
/// - **Numeric + Array**, **Array + Numeric**, and **Array + Array**:
///   broadcasted element-wise addition.
/// - **Expression + Any** and **Any + Expression**: symbolic addition.
///
/// Signature
/// ---------
/// ```zig
/// fn add(x: X, y: Y, ctx: anytype) !Add(X, Y)
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
/// `Add(@TypeOf(x), @TypeOf(y))`:
/// The result of the addition.
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
pub inline fn add(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Add(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
        .array => switch (comptime types.domain(Y)) {
            .expression => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            .array, .numeric => { // array + array, numeric + array, array + numeric
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int, .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return array.add(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .matrix => switch (comptime types.domain(Y)) {
            .expression => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            .matrix => { // matrix + matrix
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int, .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .matrix_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return matrix.add(
                    ctx.matrix_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"matrix_allocator"}),
                );
            },
            else => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .vector => switch (comptime types.domain(Y)) {
            .expression => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            .vector => { // vector + vector
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int, .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .vector_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return vector.add(
                    ctx.vector_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"vector_allocator"}),
                );
            },
            else => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domain(Y)) {
            .expression => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            .array => { // numeric + array
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int, .float, .cfloat => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                    .integer, .rational, .real, .complex => {
                        types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                                .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );
                    },
                };

                return array.add(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // numeric + numeric
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.add(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.add(x, y);
                    },
                    .cfloat => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return cfloat.add(x, y);
                    },
                    .integer => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return integer.add(ctx.allocator, x, y);
                    },
                    .rational => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return rational.add(ctx.allocator, x, y);
                    },
                    .real => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                    .complex => {
                        comptime types.validateContext(
                            @TypeOf(ctx),
                            .{
                                .allocator = .{ .type = std.mem.Allocator, .required = true },
                            },
                        );

                        return complex.add(ctx.allocator, x, y);
                    },
                }
            },
            else => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
    }
}
