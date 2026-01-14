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

/// The return type of the `log2` routine for an input of type `X`.
pub fn Log2(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Log2(types.Numeric(X))),
        .matrix => @compileError("zml.Log2 not implemented for matrices yet"),
        .vector => @compileError("zml.Log2 not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the base-2 logarithm of `x`.
///
/// The `log2` routine computes the base-2 logarithm of its input `x`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar base-2 logarithm.
/// - **Matrix**: matrix base-2 logarithm (not implemented yet).
/// - **Array**: element-wise base-2 logarithm.
/// - **Expression**: symbolic base-2 logarithm.
///
/// Signature
/// ---------
/// ```zig
/// fn log2(x: X, ctx: anytype) !Log2(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the base-2 logarithm of.
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
/// `Log2(@TypeOf(x))`:
/// The base-2 logarithm of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn log2(
    x: anytype,
    ctx: anytype,
) !Log2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.log2 not implemented for expressions yet"),
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
                else => @compileError("zml.log2 for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.log2(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.log2 not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.log2(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.log2(x);
            },
            else => @compileError("zml.log2 for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
