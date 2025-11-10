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

/// The return type of the `cbrt` routine for an input of type `X`.
pub fn Cbrt(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Cbrt(types.Numeric(X))),
        .matrix => @compileError("zml.Cbrt not implemented for matrices yet"),
        .vector => @compileError("zml.Cbrt not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the cube root of `x`, `∛x`.
///
/// The `cbrt` routine computes the cube root of its input `x`, `∛x`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar cube root.
/// - **Matrix**: matrix cube root (not implemented yet).
/// - **Array**: element-wise cube root.
/// - **Expression**: symbolic cube root.
///
/// Signature
/// ---------
/// ```zig
/// fn cbrt(x: X, ctx: anytype) !Cbrt(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the cube root of.
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
/// `Cbrt(@TypeOf(x))`:
/// The cube root of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn cbrt(
    x: anytype,
    ctx: anytype,
) !Cbrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.cbrt not implemented for expressions yet"),
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
                else => @compileError("zml.cbrt for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.cbrt(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cbrt not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.cbrt(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.cbrt(x);
            },
            else => @compileError("zml.cbrt for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
