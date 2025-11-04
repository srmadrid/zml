const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

/// The return type of the `cosh` routine for an input of type `X`.
pub fn Cosh(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Cosh(types.Numeric(X))),
        .matrix => @compileError("zml.Cosh not implemented for matrices yet"),
        .vector => @compileError("zml.Cosh not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the hyperbolic cosine of `x`.
///
/// The `cosh` routine computes the hyperbolic cosine its input `x`, validating
/// the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar hyperbolic cosine.
/// - **Matrix**: matrix hyperbolic cosine (not implemented yet).
/// - **Array**: element-wise hyperbolic cosine.
///
/// Signature
/// ---------
/// ```zig
/// fn cosh(x: X, ctx: anytype) !Cosh(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the hyperbolic cosine of.
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
/// `Cosh(@TypeOf(x))`:
/// The hyperbolic cosine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn cosh(
    x: anytype,
    ctx: anytype,
) !Cosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.cosh not defined for " ++ @typeName(X));

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
                else => @compileError("zml.cosh for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.cosh(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cosh not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.cosh(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.cosh(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.cosh(x);
            },
            else => @compileError("zml.cosh for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
