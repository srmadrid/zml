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

/// The return type of the `acos` routine for an input of type `X`.
pub fn Acos(X: type) type {
    return switch (comptime types.domainType(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Acos(types.Numeric(X))),
        .matrix => @compileError("zml.Acos not implemented for matrices yet"),
        .vector => @compileError("zml.Acos not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the arccosine of `x`, `cos⁻¹(x)`.
///
/// The `acos` routine computes the arccosine its input `x`, `cos⁻¹(x)`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar arccosine.
/// - **Matrix**: matrix arccosine (not implemented yet).
/// - **Array**: element-wise arccosine.
/// - **Expression**: symbolic arccosine.
///
/// Signature
/// ---------
/// ```zig
/// fn acos(x: X, ctx: anytype) !Acos(X)
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
/// `Acos(@TypeOf(x))`:
/// The arccosine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn acos(
    x: anytype,
    ctx: anytype,
) !Acos(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domainType(X)) {
        .expression => @compileError("zml.acos not defined for " ++ @typeName(X)),
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
                else => @compileError("zml.acos for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.acos(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.acos not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.acos(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.acos(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.acos(x);
            },
            else => @compileError("zml.acos for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
