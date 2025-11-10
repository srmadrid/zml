const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");
const expression = @import("../expression.zig");

/// The return type of the `tan` routine for an input of type `X`.
pub fn Tan(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Tan(types.Numeric(X))),
        .matrix => @compileError("zml.Tan not implemented for matrices yet"),
        .vector => @compileError("zml.Tan not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the tangent of `x`.
///
/// The `tan` routine computes the tangent its input `x`, validating the
/// provided context. It supports both fixed-precision and arbitrary-precision
/// arithmetic, as well as structured data domains. The supported domains are:
/// - **Numeric**: scalar tangent.
/// - **Matrix**: matrix tangent (not implemented yet).
/// - **Array**: element-wise tangent.
/// - **Expression**: symbolic tangent.
///
/// Signature
/// ---------
/// ```zig
/// fn tan(x: X, ctx: anytype) !Tan(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the sine of.
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
/// `Tan(@TypeOf(x))`:
/// The tangent of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn tan(
    x: anytype,
    ctx: anytype,
) !Tan(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.tan not implemented for expressions yet"),
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
                else => @compileError("zml.tan for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.tan(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.tan not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.tan(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.tan(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.tan(x);
            },
            else => @compileError("zml.tan for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
