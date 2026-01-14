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

/// The return type of the `acosh` routine for an input of type `X`.
pub fn Acosh(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Acosh(types.Numeric(X))),
        .matrix => @compileError("zml.Acosh not implemented for matrices yet"),
        .vector => @compileError("zml.Acosh not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the hyperbolic arccosine of `x`, `cosh⁻¹(x)`.
///
/// The `acosh` routine computes the hyperbolic arccosine its input `x`,
/// `cosh⁻¹(x)`, validating the provided context. It supports both
/// fixed-precision and arbitrary-precision arithmetic, as well as structured
/// data domains. The supported domains are:
/// - **Numeric**: scalar hyperbolic arccosine.
/// - **Matrix**: matrix hyperbolic arccosine (not implemented yet).
/// - **Array**: element-wise hyperbolic arccosine.
/// - **Expression**: symbolic hyperbolic arccosine.
///
/// Signature
/// ---------
/// ```zig
/// fn acosh(x: X, ctx: anytype) !Acosh(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the hyperbolic arccosine of.
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
/// `Acosh(@TypeOf(x))`:
/// The hyperbolic arccosine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn acosh(
    x: anytype,
    ctx: anytype,
) !Acosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.acosh for " ++ @typeName(X) ++ " not implemented yet"),
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
                else => @compileError("zml.acosh for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.acosh(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.acosh not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.acosh(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.acosh(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.acosh(x);
            },
            else => @compileError("zml.acosh for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
