const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");
const expression = @import("../expression.zig");

/// The return type of the `sqrt` routine for an input of type `X`.
pub fn Sqrt(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Sqrt(types.Numeric(X))),
        .matrix => @compileError("zml.Sqrt not implemented for matrices yet"),
        .vector => @compileError("zml.Sqrt not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the square root of `x`, `√x`.
///
/// The `sqrt` routine computes the square root of its input `x`, `√x`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar square root.
/// - **Matrix**: matrix square root (not implemented yet).
/// - **Array**: element-wise square root.
/// - **Expression**: symbolic square root.
///
/// Signature
/// ---------
/// ```zig
/// fn sqrt(x: X, ctx: anytype) !Sqrt(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the square root of.
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
/// `Sqrt(@TypeOf(x))`:
/// The square root of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn sqrt(
    x: anytype,
    ctx: anytype,
) !Sqrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.sqrt not implemented for expressions yet"),
        .array => {
            comptime switch (types.numericType(types.Numeric(X))) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
                else => @compileError("zml.sqrt for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.sqrt(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sqrt not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sqrt(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sqrt(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sqrt(x);
            },
            else => @compileError("zml.sqrt for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
