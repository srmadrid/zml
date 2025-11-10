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

/// The return type of the `exp` routine for an input of type `X`.
pub fn Exp(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Exp(types.Numeric(X))),
        .matrix => @compileError("zml.Exp not implemented for matrices yet"),
        .vector => @compileError("zml.Exp not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the exponential of `x`, `eˣ`.
///
/// The `exp` routine computes the exponential of its input `x`, `eˣ`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar exponential.
/// - **Matrix**: matrix exponential (not implemented yet).
/// - **Array**: element-wise exponential.
/// - **Expression**: symbolic exponential.
///
/// Signature
/// ---------
/// ```zig
/// fn exp(x: X, ctx: anytype) !Exp(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the exponential of.
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
/// `Exp(@TypeOf(x))`:
/// The exponential of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn exp(
    x: anytype,
    ctx: anytype,
) !Exp(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.exp not implemented for " ++ @typeName(X) ++ " yet"),
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
                else => @compileError("zml.exp for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.exp(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .matrix => @compileError("zml.exp not implemented for " ++ @typeName(X) ++ " yet"),
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.exp(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.exp(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.exp(x);
            },
            else => @compileError("zml.exp for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
