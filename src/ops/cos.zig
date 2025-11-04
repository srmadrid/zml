const std = @import("std");

const types = @import("../types.zig");
const EnsureArray = types.EnsureArray;
const EnsureFloat = types.EnsureFloat;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");

const array = @import("../array.zig");

/// The return type of the `cos` routine for an input of type `X`.
pub fn Cos(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Cos(types.Numeric(X))),
        .matrix => @compileError("zml.Cos not implemented for matrices yet"),
        .vector => @compileError("zml.Cos not defined for " ++ @typeName(X)),
        .numeric => types.EnsureFloat(X),
    };
}

/// Returns the cosine of `x`.
///
/// The `cos` routine computes the cosine its input `x`, validating the provided
/// context. It supports both fixed-precision and arbitrary-precision
/// arithmetic, as well as structured data domains. The supported domains are:
/// - **Numeric**: scalar cosine.
/// - **Matrix**: matrix cosine (not implemented yet).
/// - **Array**: element-wise cosine.
///
/// Signature
/// ---------
/// ```zig
/// fn cos(x: X, ctx: anytype) !Cos(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the cosine of.
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
/// `Cos(@TypeOf(x))`:
/// The cosine of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision or a structured data type.
pub inline fn cos(
    x: anytype,
    ctx: anytype,
) !Cos(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isNumeric(X))
        @compileError("zml.cos not defined for " ++ @typeName(X));

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
                else => @compileError("zml.cos for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.cos(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.cos not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.cos(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.cos(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.cos(x);
            },
            else => @compileError("zml.cos for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
