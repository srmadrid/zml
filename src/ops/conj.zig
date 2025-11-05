const std = @import("std");

const types = @import("../types.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// The return type of the `conj` routine for an input of type `X`.
pub fn Conj(X: type) type {
    return switch (comptime types.domainType(X)) {
        .array => types.EnsureArray(X, Conj(types.Numeric(X))),
        .matrix => @compileError("zml.Conj not implemented for matrices yet"),
        .vector => @compileError("zml.Conj not implemented for vectors yet"),
        .numeric => X,
    };
}

/// Returns the complex conjugate of `x`.
///
/// The `conj` routine computes the complex conjugate of its input `x`,
/// validating the provided context. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported domains are:
/// - **Numeric**: scalar complex conjugate.
/// - **Vector**: element-wise complex conjugate.
/// - **Matrix**: element-wise complex conjugate.
/// - **Array**: element-wise complex conjugate.
///
/// Signature
/// ---------
/// ```zig
/// fn conj(x: X, ctx: anytype) !Conj(X)
/// ```
///
/// Parameters
/// ----------
/// `x` (`anytype`):
/// The operand to compute the complex conjugate of.
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
/// `Conj(@TypeOf(x))`:
/// The complex conjugate of `x`.
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
pub inline fn conj(
    x: anytype,
    ctx: anytype,
) !Conj(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isMatrix(X) and
        !types.isVector(X) and !types.isNumeric(X))
        @compileError("zml.conj not defined for " ++ @typeName(X));

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
                else => @compileError("zml.conj for " ++ @typeName(X) ++ " not implemented yet"),
            };

            return array.conj(
                ctx.array_allocator,
                x,
                types.stripStruct(ctx, &.{"array_allocator"}),
            );
        },
        .numeric => switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.conj not defined for " ++ @typeName(X)),
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

                return x.conj();
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
            .real => @compileError("zml.conj not implemented for " ++ @typeName(X) ++ " yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false },
                    },
                );

                return try x.conj(types.getFieldOrDefault(ctx, "allocator", ?std.mem.Allocator, null));
            },
            .expression => @compileError("zml.conj not implemented for " ++ @typeName(X) ++ " yet"),
        },
        else => unreachable,
    }
}
