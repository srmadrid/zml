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

/// The return type of the `erfc` routine for an input of type `X`.
pub fn Erfc(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Erfc(types.Numeric(X))),
        .matrix => @compileError("zml.Erfc not implemented for matrices yet"),
        .vector => @compileError("zml.Erfc not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the complementary error function at `x`, `1 - ²/π ∫₀ˣ e^(−t²)dt`.
///
/// The `erfc` routine computes the complementary error function at its input
/// `x`, `1 - ²/π ∫₀ˣ e^(−t²)dt`, validating the provided context. It supports
/// both fixed-precision and arbitrary-precision arithmetic, as well as
/// structured data domains. The supported domains are:
/// - **Numeric**: scalar complementary error function.
/// - **Matrix**: matrix complementary error function (not implemented yet).
/// - **Array**: element-wise complementary error function.
/// - **Expression**: symbolic complementary error function.
///
/// Signature
/// ---------
/// ```zig
/// fn erfc(x: X, ctx: anytype) !Erfc(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the complementary error function at.
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
/// `Erfc(@TypeOf(x))`:
/// The complementary error function at `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn erfc(
    x: anytype,
    ctx: anytype,
) !Erfc(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.erfc not implemented for expressions yet"),
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
                else => @compileError("zml.erfc for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.erfc(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.erfc not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.erfc(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.erfc(x);
            },
            else => @compileError("zml.erfc for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
