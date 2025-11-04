const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

/// The return type of the `exp10` routine for an input of type `X`.
pub fn Exp10(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Exp10(types.Numeric(X))),
        .matrix => @compileError("zml.Exp10 not implemented for matrices yet"),
        .vector => @compileError("zml.Exp10 not defined for " ++ @typeName(X)),
        .numeric => X,
    };
}

/// Returns the base-10 exponential of `x`, `10ˣ`.
///
/// The `exp10` routine computes the base-10 exponential of its input `x`,
/// `10ˣ`, validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar base-10 exponential.
/// - **Matrix**: matrix base-10 exponential (not implemented yet).
/// - **Array**: element-wise base-10 exponential.
///
/// Signature
/// ---------
/// ```zig
/// fn exp10(x: X, ctx: anytype) !Exp10(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the base-10 exponential of.
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
/// `Exp10(@TypeOf(x))`:
/// The base-10 exponential of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn exp10(
    x: anytype,
    ctx: anytype,
) !Exp10(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.exp10 not defined for " ++ @typeName(X));

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
                else => @compileError("zml.exp10 for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.exp10(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.exp10 not defined for " ++ @typeName(X)),
            .int => @compileError("zml.exp10 not implemented for " ++ @typeName(X) ++ " yet"),
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.exp10(x);
            },
            else => @compileError("zml.exp10 for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
