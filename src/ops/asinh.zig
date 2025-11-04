const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

/// The return type of the `asinh` routine for an input of type `X`.
pub fn Asinh(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Asinh(types.Numeric(X))),
        .matrix => @compileError("zml.Asinh not implemented for matrices yet"),
        .vector => @compileError("zml.Asinh not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the hyperbolic arcsine of `x`, `sinh⁻¹(x)`.
///
/// The `asinh` routine computes the hyperbolic arcsine its input `x`,
/// `sinh⁻¹(x)`, validating the provided context. It supports both
/// fixed-precision and arbitrary-precision arithmetic, as well as structured
/// data domains. The supported domains are:
/// - **Numeric**: scalar hyperbolic arcsine.
/// - **Matrix**: matrix hyperbolic arcsine (not implemented yet).
/// - **Array**: element-wise hyperbolic arcsine.
///
/// Signature
/// ---------
/// ```zig
/// fn asinh(x: X, ctx: anytype) !Asinh(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the hyperbolic arcsine of.
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
/// `Asinh(@TypeOf(x))`:
/// The hyperbolic arcsine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn asinh(
    x: anytype,
    ctx: anytype,
) !Asinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.asinh not defined for " ++ @typeName(X));

    switch (comptime types.domainType(X)) {
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
                else => @compileError("zml.asinh for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.asinh(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.asinh not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.asinh(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.asinh(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.asinh(x);
            },
            else => @compileError("zml.asinh for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
