const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

/// The return type of the `atan2` routine for inputs of types `X` and `Y`.
pub fn Atan2(X: type, Y: type) type {
    return switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array => types.EnsureArray(Y, Atan2(types.Numeric(X), types.Numeric(Y))),
            .matrix => @compileError("zml.Atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .vector => @compileError("zml.Atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.EnsureArray(Y, Atan2(types.Numeric(X), Y)),
        },
        .matrix => @compileError("zml.Atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .vector => @compileError("zml.Atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .numeric => switch (comptime types.domainType(Y)) {
            .array => types.EnsureArray(Y, Atan2(X, types.Numeric(Y))),
            .matrix => @compileError("zml.Atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .vector => @compileError("zml.Atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
            .numeric => types.EnsureFloat(types.Coerce(X, Y)),
        },
    };
}

/// Performs the arctangent of `y/x` using the signs of both arguments to
/// determine the correct quadrant.
///
/// The `atan2` routine computes the arctangent of `y/x`, using the signs of
/// both arguments to determine the correct quadrant, automatically coercing
/// compatible operand types and validating the provided context. The operation
/// is performed in the coerced precision of the operands, and the resulting
/// value is returned as a new value. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric ^ Numeric**: scalar arctangent.
/// - **Numeric * Array**, **Array ^ Numeric**, and **Array ^ Array**:
///   broadcasted element-wise arctangent.
///
/// Signature
/// ---------
/// ```zig
/// fn atan2(x: X, y: Y, ctx: anytype) !Atan2(X, Y)
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
/// `Atan2(@TypeOf(x), @TypeOf(y))`:
/// The result of the arctangent.
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
pub inline fn atan2(
    x: anytype,
    y: anytype,
    ctx: anytype,
) !Atan2(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y) and
        !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    const C: type = types.Coerce(X, Y);

    switch (comptime types.domainType(X)) {
        .array => switch (comptime types.domainType(Y)) {
            .array, .numeric => { // atan2(array, array), atan2(array, numeric)
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

                return array.atan2(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => @compileError("zml.atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        .numeric => switch (comptime types.domainType(Y)) {
            .array => { // atan2(numeric, array)
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

                return array.atan2(
                    ctx.array_allocator,
                    x,
                    y,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            .numeric => { // atan2(numeric, numeric)
                switch (comptime types.numericType(C)) {
                    .bool => @compileError("zml.atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                    .int => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.atan2(x, y);
                    },
                    .float => {
                        comptime types.validateContext(@TypeOf(ctx), .{});

                        return float.atan2(x, y);
                    },
                    else => @compileError("zml.atan2 between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
                }
            },
            else => @compileError("zml.atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        },
        else => @compileError("zml.atan2 not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
    }
}
