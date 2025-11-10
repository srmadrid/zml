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

/// The return type of the `sin` routine for an input of type `X`.
pub fn Sin(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Sin(types.Numeric(X))),
        .matrix => @compileError("zml.Sin not implemented for matrices yet"),
        .vector => @compileError("zml.sin not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the sine of `x`.
///
/// The `sin` routine computes the sine its input `x`, validating the provided
/// context. It supports both fixed-precision and arbitrary-precision
/// arithmetic, as well as structured data domains. The supported domains are:
/// - **Numeric**: scalar sine.
/// - **Matrix**: matrix sine (not implemented yet).
/// - **Array**: element-wise sine.
/// - **Expression**: symbolic sine.
///
/// Signature
/// ---------
/// ```zig
/// fn sin(x: X, ctx: anytype) !Sin(X)
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
/// `Sin(@TypeOf(x))`:
/// The sine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn sin(
    x: anytype,
    ctx: anytype,
) !Sin(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.sin not implemented for expressions yet"),
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
                else => @compileError("zml.sin for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.sin(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.sin not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sin(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sin(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sin(x);
            },
            else => @compileError("zml.sin for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
