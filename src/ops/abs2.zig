const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// The return type of the `abs2` routine for an input of type `X`.
pub fn Abs2(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Abs2(types.Numeric(X))),
        .matrix => @compileError("zml.Abs2 not implemented for matrices yet"),
        .vector => types.Scalar(types.Numeric(X)),
        .numeric => types.Scalar(X),
    };
}

/// Returns the squared absolute value of `x`.
///
/// The `abs2` routine computes the absolute value of its input `x`, validating
/// the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: squared scalar absolute value.
/// - **Vector**: not defined yet, eventually |x₁|² + |x₂|² + ... + |xₙ|².
/// - **Matrix**: not defined yet, eventually AᴴA.
/// - **Array**: element-wise squared absolute value.
///
/// Signature
/// ---------
/// ```zig
/// fn abs2(x: X, ctx: anytype) !Abs2(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the squared absolute value of.
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
/// `Abs2(@TypeOf(x))`:
/// The squared absolute value of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the type is of arbitrary
/// precision and an allocator is provided, or a structured data type.
///
/// Notes
/// -----
/// For some arbitrary-precision numeric types, providing an allocator in the
/// context is optional. If provided, a new value will be allocated for the
/// result. If not provided, the operation will return a view.
pub inline fn abs2(
    x: anytype,
    ctx: anytype,
) !Abs2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isMatrix(X) and
        !types.isVector(X) and !types.isNumeric(X))
        @compileError("zml.abs2 not defined for " ++ @typeName(X));

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
                .integer, .rational => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
                        },
                    );
                },
                .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
                else => @compileError("zml.abs2 for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.abs2(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .matrix => @compileError("zml.abs2 for " ++ @typeName(X) ++ " not implemented yet"),
        .vector => @compileError("zml.abs2 for " ++ @typeName(X) ++ " not implemented yet"),
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs2 not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.mul(x, x, .default);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.mul(x, x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.abs2(x);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return integer.mul(ctx.allocator, x, x);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return rational.mul(ctx.allocator, x, x);
            },
            .real => @compileError("zml.abs2 for " ++ @typeName(X) ++ " not implemented yet"),
            .complex => @compileError("zml.abs2 for " ++ @typeName(X) ++ " not implemented yet"),
            .expression => @compileError("zml.abs2 for " ++ @typeName(X) ++ " not implemented yet"),
        },
    }
}
