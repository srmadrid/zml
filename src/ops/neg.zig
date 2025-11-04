const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const array = @import("../array.zig");

/// The return type of the `neg` routine for an input of type `X`.
pub fn Neg(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, types.Numeric(X)),
        .matrix => X,
        .vector => X,
        .numeric => X,
    };
}

/// Returns `x` negated.
///
/// The `neg` routine computes the negation of its input `x`, validating the
/// provided context. It supports both fixed-precision and arbitrary-precision
/// arithmetic, as well as structured data domains. The supported domains are:
/// - **Numeric**: scalar negation.
/// - **Vector**: element-wise negation.
/// - **Matrix**: element-wise negation.
/// - **Array**: element-wise negation.
///
/// Signature
/// ---------
/// ```zig
/// fn neg(x: X, ctx: anytype) !Neg(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the negation of.
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
/// `Neg(@TypeOf(x))`:
/// The negation of `x`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`:
/// If memory allocation fails. Can only happen if the coerced type is of
/// arbitrary precision and an allocator is provided, or a structured data type.
///
/// Notes
/// -----
/// For some arbitrary-precision numeric types, providing an allocator in the
/// context is optional. If provided, a new value will be allocated for the
/// result. If not provided, the operation will return a view.
pub inline fn neg(
    x: anytype,
    ctx: anytype,
) !Neg(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isMatrix(X) and
        !types.isVector(X) and !types.isNumeric(X))
        @compileError("zml.neg not defined for " ++ @typeName(X));

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
                .integer, .rational, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .element_allocator = .{ .type = ?std.mem.Allocator, .required = false },
                        },
                    );
                },
                else => @compileError("zml.neg for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.neg(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .matrix => @compileError("zml.neg for " ++ @typeName(X) ++ " not implemented yet"),
        .vector => @compileError("zml.neg for " ++ @typeName(X) ++ " not implemented yet"),
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.neg not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return -x;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return -x;
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.neg(x);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return integer.neg(types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null), x);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return rational.neg(types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null), x);
            },
            .real => @compileError("zml.neg for " ++ @typeName(X) ++ " not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return complex.neg(ctx.allocator, x);
            },
            .expression => @compileError("zml.neg for " ++ @typeName(X) ++ " not implemented yet"),
        },
    }
}
