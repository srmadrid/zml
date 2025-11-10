const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");
const expression = @import("../expression.zig");

/// The return type of the `tanh` routine for an input of type `X`.
pub fn Tanh(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Tanh(types.Numeric(X))),
        .matrix => @compileError("zml.Tanh not implemented for matrices yet"),
        .vector => @compileError("zml.Tanh not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the hyperbolic tangent of `x`.
///
/// The `tanh` routine computes the hyperbolic tangent its input `x`, validating
/// the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar hyperbolic tangent.
/// - **Matrix**: matrix hyperbolic tangent (not implemented yet).
/// - **Array**: element-wise hyperbolic tangent.
/// - **Expression**: symbolic hyperbolic tangent.
///
/// Signature
/// ---------
/// ```zig
/// fn tanh(x: X, ctx: anytype) !Tanh(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the hyperbolic tangent of.
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
/// `Tanh(@TypeOf(x))`:
/// The hyperbolic tangent of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn tanh(
    x: anytype,
    ctx: anytype,
) !Tanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.tanh not implemented for expressions yet"),
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
                else => @compileError("zml.tanh for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.tanh(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.tanh not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.tanh(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.tanh(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.tanh(x);
            },
            else => @compileError("zml.tanh for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
