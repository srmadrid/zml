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

/// The return type of the `pow` routine for inputs of types `X` and `Y`.
pub fn Pow(X: type, Y: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => switch (comptime types.domainType(Y)) {
            .array => expression.Expression,
            .matrix => @compileError("zml.Pow not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet"),
            .vector => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.EnsureArray(Y, Pow(types.Numeric(X), Y)),
        },
        .array => switch (comptime types.domainType(Y)) {
            .expression => expression.Expression,
            .array => types.EnsureArray(Y, Pow(types.Numeric(X), types.Numeric(Y))),
            .matrix => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .vector => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.EnsureArray(Y, Pow(types.Numeric(X), Y)),
        },
        .matrix => switch (comptime types.domainType(Y)) {
            .expression => expression.Expression,
            .array => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .matrix => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .vector => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => @compileError("zml.Pow not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet"),
        },
        .vector => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .numeric => switch (comptime types.domainType(Y)) {
            .expression => expression.Expression,
            .array => types.EnsureArray(Y, Pow(X, types.Numeric(Y))),
            .matrix => @compileError("zml.Pow not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet"),
            .vector => @compileError("zml.Pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.Coerce(X, Y),
        },
    };
}

/// Performs the exponentiation operation `xʸ`.
///
/// The `pow` routine computes the exponentiation `xʸ`, automatically coercing
/// compatible operand types and validating the provided context. The operation
/// is performed in the coerced precision of the operands, and the resulting
/// value is returned as a new value. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric ^ Numeric**: scalar exponentiation.
/// - **Numeric ^ Matrix** and **Matrix ^ Numeric**: matrix exponentiation.
/// - **Numeric * Array**, **Array ^ Numeric**, and **Array ^ Array**:
///   broadcasted element-wise exponentiation.
/// - **Numeric ^ Expression**, **Matrix ^ Expression**, **Array ^ Expression**,
///   **Expression ^ Numeric**, **Expression ^ Matrix**, **Expression ^ Array**,
///   and **Expression ^ Expression**: symbolic exponentiation.
///
/// Signature
/// ---------
/// ```zig
/// fn pow(x: X, y: Y, ctx: anytype) !Pow(X, Y)
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
/// `Pow(@TypeOf(x), @TypeOf(y))`:
/// The result of the exponentiation.
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
///
/// `int.Error.NegativeExponent`:
/// If both operands are integers and the exponent is negative.
pub inline fn pow(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Pow(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .expression => @compileError("zml.pow not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet"),
            .array, .numeric => { // array^array, array^numeric
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int, .float, .cfloat => {
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

                return array.pow(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .expression => @compileError("zml.pow not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet"),
            .array => { // numeric^array
                comptime switch (types.numericType(types.Numeric(C))) {
                    .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int, .float, .cfloat => {
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

                return array.pow(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // numeric^numeric
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return int.pow(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.pow(x, y);
                    },
                    .cfloat => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return cfloat.pow(x, y);
                    },
                    else => @compileError("zml.pow between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.pow not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
