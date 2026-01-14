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

/// The return type of the `atanh` routine for an input of type `X`.
pub fn Atanh(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Atanh(types.Numeric(X))),
        .matrix => @compileError("zml.Atanh not implemented for matrices yet"),
        .vector => @compileError("zml.Atanh not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the hyperbolic arctangent of `x`, `tanh⁻¹(x)`.
///
/// The `atanh` routine computes the hyperbolic arctangent of its input `x`,
/// `tanh⁻¹(x)`, validating the provided context. It supports both
/// fixed-precision and arbitrary-precision arithmetic, as well as structured
/// data domains. The supported domains are:
/// - **Numeric**: scalar hyperbolic arctangent.
/// - **Matrix**: matrix hyperbolic arctangent (not implemented yet).
/// - **Array**: element-wise hyperbolic arctangent.
/// - **Expression**: symbolic hyperbolic arctangent.
///
/// Signature
/// ---------
/// ```zig
/// fn atanh(x: X, ctx: anytype) !Atanh(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the hyperbolic arctangent of.
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
/// `Atanh(@TypeOf(x))`:
/// The hyperbolic arctangent of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn atanh(
    x: anytype,
    ctx: anytype,
) !Atanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.atanh not implemented for expressions yet"),
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
                else => @compileError("zml.atanh for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.atanh(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.atanh not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.atanh(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.atanh(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.atanh(x);
            },
            else => @compileError("zml.atanh for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
