const std = @import("std");

const types = @import("../types.zig");
const Numeric = types.Numeric;
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// The return type of the `re` routine for an input of type `X`.
pub fn Re(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Re(types.Numeric(X))),
        .matrix => @compileError("zml.Re not implemented for matrices yet"),
        .vector => @compileError("zml.Re not implemented for vectors yet"),
        .numeric => types.Scalar(X),
    };
}

/// Returns the real part of `x`.
///
/// The `re` routine computes the real part of its input `x`, validating the
/// provided context. It supports both fixed-precision and arbitrary-precision
/// arithmetic, as well as structured data domains. The supported domains are:
/// - **Numeric**: scalar real part.
/// - **Vector**: element-wise real part.
/// - **Matrix**: element-wise real part.
/// - **Array**: element-wise real part.
///
/// Signature
/// ---------
/// ```zig
/// fn re(x: X, ctx: anytype) !Re(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the real part of.
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
/// `Re(@TypeOf(x))`:
/// The real part of `x`.
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
pub inline fn re(
    x: anytype,
    ctx: anytype,
) !Re(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isMatrix(X) and
        !types.isVector(X) and !types.isNumeric(X))
        @compileError("zml.re not defined for " ++ @typeName(X));

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
                else => @compileError("zml.re for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.re(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .matrix => @compileError("zml.re not implemented for matrices yet"),
        .vector => @compileError("zml.re not implemented for vectors yet"),
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.re not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return x;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return x;
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return x.re;
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                if (types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null)) |allocator| {
                    return x.copy(allocator);
                } else {
                    var r: integer.Integer = x;
                    r.flags.owns_data = false;
                    return r;
                }
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                if (types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null)) |allocator| {
                    return x.copy(allocator);
                } else {
                    var r: rational.Rational = x;
                    r.flags.owns_data = false;
                    return r;
                }
            },
            .real => @compileError("zml.re for real numbers not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                if (types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null)) |allocator| {
                    return x.re.copy(allocator);
                } else {
                    var r: types.Scalar(X) = x.re;
                    r.flags.owns_data = false;
                    return r;
                }
            },
            .expression => @compileError("zml.re for expression types not implemented yet"),
        },
    }
}
