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

/// The return type of the `erf` routine for an input of type `X`.
pub fn Erf(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Erf(types.Numeric(X))),
        .matrix => @compileError("zml.Erf not implemented for matrices yet"),
        .vector => @compileError("zml.Erf not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the error function at `x`, `²/π ∫₀ˣ e^(−t²)dt`.
///
/// The `erf` routine computes the error function at its input `x`,
/// `²/π ∫₀ˣ e^(−t²)dt`, validating the provided context. It supports both
/// fixed-precision and arbitrary-precision arithmetic, as well as structured
/// data domains. The supported domains are:
/// - **Numeric**: scalar error function.
/// - **Matrix**: matrix error function (not implemented yet).
/// - **Array**: element-wise error function.
/// - **Expression**: symbolic error function.
///
/// Signature
/// ---------
/// ```zig
/// fn erf(x: X, ctx: anytype) !Erf(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the error function at.
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
/// `Erf(@TypeOf(x))`:
/// The error function at `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn erf(
    x: anytype,
    ctx: anytype,
) !Erf(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.erf not implemented for expressions yet"),
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
                else => @compileError("zml.erf for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.erf(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.erf not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.erf(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.erf(x);
            },
            else => @compileError("zml.erf for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
