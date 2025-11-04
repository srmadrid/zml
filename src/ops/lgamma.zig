const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

/// The return type of the `lgamma` routine for an input of type `X`.
pub fn Lgamma(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Lgamma(types.Numeric(X))),
        .matrix => @compileError("zml.Lgamma not implemented for matrices yet"),
        .vector => @compileError("zml.Lgamma not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the log-gamma function at `x`, `log(∫₀^∞ tˣ⁻¹ e⁻ᵗ dt)`.
///
/// The `lgamma` routine computes the log-gamma function at its input `x`,
/// `log(∫₀^∞ tˣ⁻¹ e⁻ᵗ dt)`, validating the provided context. It supports both
/// fixed-precision and arbitrary-precision arithmetic, as well as structured
/// data domains. The supported domains are:
/// - **Numeric**: scalar log-gamma function.
/// - **Matrix**: matrix log-gamma function (not implemented yet).
/// - **Array**: element-wise log-gamma function.
///
/// Signature
/// ---------
/// ```zig
/// fn lgamma(x: X, ctx: anytype) !Lgamma(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the log-gamma function at.
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
/// `Lgamma(@TypeOf(x))`:
/// The log-gamma function at `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn lgamma(
    x: anytype,
    ctx: anytype,
) !Lgamma(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.lgamma not defined for " ++ @typeName(X));

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
                else => @compileError("zml.lgamma for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.lgamma(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.lgamma not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.lgamma(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.lgamma(x);
            },
            else => @compileError("zml.lgamma for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
