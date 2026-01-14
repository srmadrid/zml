const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");
const expression = @import("../expression.zig");

/// The return type of the `abs1` routine for an input of type `X`.
pub fn Abs1(X: type) type {
    return switch (comptime types.domain(X)) {
        .expression => expression.Expression,
        .array => types.EnsureArray(X, Abs1(types.Numeric(X))),
        .matrix => @compileError("zml.Abs1 not defined for " ++ @typeName(X)),
        .vector => @compileError("zml.Abs1 not defined for " ++ @typeName(X)),
        .numeric => types.Scalar(X),
    };
}

/// Returns the absolute value of `x` for real inputs, or the sum of the
/// absolute values of the real and imaginary parts for complex inputs.
///
/// The `abs1` routine computes the absolute value of its input `x`, validating
/// the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar absolute value.
/// - **Array**: element-wise absolute value.
/// - **Expression**: symbolic absolute value.
///
/// Signature
/// ---------
/// ```zig
/// fn abs1(x: X, ctx: anytype) !Abs1(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the absolute value of.
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
/// `Abs1(@TypeOf(x))`:
/// The absolute value of `x` for real inputs, or the sum of the absolute values
/// of the real and imaginary parts for complex inputs.
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
pub inline fn abs1(
    x: anytype,
    ctx: anytype,
) !Abs1(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (comptime types.domain(X)) {
        .expression => @compileError("zml.abs1 for " ++ @typeName(X) ++ " not implemented yet"),
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
                            .element_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
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
                else => @compileError("zml.abs1 for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.abs1(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.abs1 not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.abs(x);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.abs(x);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.abs1(x);
            },
            .integer => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                return integer.abs(types.getFieldOrDefault(ctx, spec, "allocator"), x);
            },
            .rational => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                return rational.abs(types.getFieldOrDefault(ctx, spec, "allocator"), x);
            },
            .real => @compileError("zml.abs1 for " ++ @typeName(X) ++ " not implemented yet"),
            .complex => @compileError("zml.abs1 for " ++ @typeName(X) ++ " not implemented yet"),
        },
        else => unreachable,
    }
}
