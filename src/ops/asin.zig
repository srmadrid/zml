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

/// The return type of the `asin` routine for an input of type `X`.
pub fn Asin(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Asin(types.Numeric(X))),
        .matrix => @compileError("zml.Asin not implemented for matrices yet"),
        .vector => @compileError("zml.Asin not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the arcsine of `x`, `sin⁻¹(x)`.
///
/// The `asin` routine computes the arcsine its input `x`, `sin⁻¹(x)`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar arcsine.
/// - **Matrix**: matrix arcsine (not implemented yet).
/// - **Array**: element-wise arcsine.
/// - **Expression**: symbolic arcsine.
///
/// Signature
/// ---------
/// ```zig
/// fn asin(x: X, ctx: anytype) !Asin(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the arcsine of.
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
/// `Asin(@TypeOf(x))`:
/// The arcsine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn asin(
    x: anytype,
    ctx: anytype,
) !Asin(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.asin not implemented for expressions yet"),
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
                else => @compileError("zml.asin for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.asin(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.asin not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.asin(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.asin(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.asin(x);
            },
            else => @compileError("zml.asin for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
