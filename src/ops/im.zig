const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// The return type of the `im` routine for an input of type `X`.
pub fn Im(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Im(types.Numeric(X))),
        .matrix => @compileError("zml.Im not implemented for matrices yet"),
        .vector => @compileError("zml.Im not implemented for vectors yet"),
        .numeric => types.Scalar(X),
    };
}

/// Returns the imaginary part of `x`.
///
/// The `im` routine computes the imaginary part of its input `x`, validating
/// the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar imaginary part.
/// - **Vector**: element-wise imaginary part.
/// - **Matrix**: element-wise imaginary part.
/// - **Array**: element-wise imaginary part.
///
/// Signature
/// ---------
/// ```zig
/// fn im(x: X, ctx: anytype) !Im(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the imaginary part of.
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
/// `Im(@TypeOf(x))`:
/// The imaginary part of `x`.
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
pub inline fn im(
    x: anytype,
    ctx: anytype,
) !Im(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isMatrix(X) and
        !types.isVector(X) and !types.isNumeric(X))
        @compileError("zml.im not defined for " ++ @typeName(X));

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
                else => @compileError("zml.im for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.im(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .matrix => @compileError("zml.im not implemented for matrices yet"),
        .vector => @compileError("zml.im not implemented for vectors yet"),
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.im not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return 0;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return 0.0;
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return x.im;
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return constants.zero(integer.Integer, ctx);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return constants.zero(rational.Rational, ctx);
            },
            .real => @compileError("zml.im for real numbers not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                if (types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null)) |allocator| {
                    return x.im.copy(allocator);
                } else {
                    var r: types.Scalar(X) = x.im;
                    r.flags.owns_data = false;
                    return r;
                }
            },
            .expression => @compileError("zml.im for expression types not implemented yet"),
        },
    }
}
