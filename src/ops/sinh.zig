const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");
const expression = @import("../expression.zig");

/// The return type of the `sinh` routine for an input of type `X`.
pub fn Sinh(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Sinh(types.Numeric(X))),
        .matrix => @compileError("zml.Sinh not implemented for matrices yet"),
        .vector => @compileError("zml.sinh not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the hyperbolic sine of `x`.
///
/// The `sinh` routine computes the hyperbolic sine its input `x`, validating
/// the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar hyperbolic sine.
/// - **Matrix**: matrix hyperbolic sine (not implemented yet).
/// - **Array**: element-wise hyperbolic sine.
/// - **Expression**: symbolic hyperbolic sine.
///
/// Signature
/// ---------
/// ```zig
/// fn sinh(x: X, ctx: anytype) !Sinh(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the hyperbolic sine of.
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
/// `Sinh(@TypeOf(x))`:
/// The hyperbolic sine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn sinh(
    x: anytype,
    ctx: anytype,
) !Sinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.sinh not implemented for expressions yet"),
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
                else => @compileError("zml.sinh for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.sinh(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sinh not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sinh(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sinh(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sinh(x);
            },
            else => @compileError("zml.sinh for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
