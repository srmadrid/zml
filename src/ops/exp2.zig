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

/// The return type of the `exp2` routine for an input of type `X`.
pub fn Exp2(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Exp2(types.Numeric(X))),
        .matrix => @compileError("zml.Exp2 not implemented for matrices yet"),
        .vector => @compileError("zml.Exp2 not defined for " ++ @typeName(X)),
        .numeric => X,
    };
}

/// Returns the base-2 exponential of `x`, `2ˣ`.
///
/// The `exp2` routine computes the base-2 exponential of its input `x`,
/// `2ˣ`, validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar base-2 exponential.
/// - **Matrix**: matrix base-2 exponential (not implemented yet).
/// - **Array**: element-wise base-2 exponential.
/// - **Expression**: symbolic base-2 exponential.
///
/// Signature
/// ---------
/// ```zig
/// fn exp2(x: X, ctx: anytype) !Exp2(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the base-2 exponential of.
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
/// `Exp2(@TypeOf(x))`:
/// The base-2 exponential of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn exp2(
    x: anytype,
    ctx: anytype,
) !Exp2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.exp2 not implemented for expressions yet"),
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
                else => @compileError("zml.exp2 for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.exp2(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp2 not defined for " ++ @typeName(X)),
            .int => @compileError("zml.exp2 not implemented for " ++ @typeName(X) ++ " yet"),
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.exp2(x);
            },
            else => @compileError("zml.exp2 for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
