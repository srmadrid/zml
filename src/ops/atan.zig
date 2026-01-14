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

/// The return type of the `atan` routine for an input of type `X`.
pub fn Atan(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Atan(types.Numeric(X))),
        .matrix => @compileError("zml.Atan not implemented for matrices yet"),
        .vector => @compileError("zml.Atan not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the arctangent of `x`, `tan⁻¹(x)`.
///
/// The `atan` routine computes the arctangent its input `x`, `tan⁻¹(x)`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar arctangent.
/// - **Matrix**: matrix arctangent (not implemented yet).
/// - **Array**: element-wise arctangent.
/// - **Expression**: symbolic arctangent.
///
/// Signature
/// ---------
/// ```zig
/// fn atan(x: X, ctx: anytype) !Atan(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the arctangent of.
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
/// `Atan(@TypeOf(x))`:
/// The arctangent of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn atan(
    x: anytype,
    ctx: anytype,
) !Atan(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.atan not implemented for expressions yet"),
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
                else => @compileError("zml.atan for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.atan(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.atan not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.atan(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.atan(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.atan(x);
            },
            else => @compileError("zml.atan for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
