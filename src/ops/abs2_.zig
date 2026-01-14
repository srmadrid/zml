const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// Performs in-place computation of the squared absolute value of `x`.
///
/// The `abs2_` routine computes the squared absolute value of its input `x`, and stores
/// the result directly into `o`, automatically validating the provided
/// context. The operation is performed in the input's precision, and the
/// result is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = abs2(Numeric)**: squared scalar absolute value.
/// - **Numeric = abs2(Vector)**: not defined yet, eventually |x₁|² + |x₂|² + ... + |xₙ|².
/// - **Matrix = abs2(Matrix)**: not defined yet, eventually AᴴA.
/// - **Array = abs2(Array)**: element-wise squared absolute value.
/// - **Expression = abs2(Expression)**: symbolic squared absolute value.
///
/// Signature
/// ---------
/// ```zig
/// fn abs2_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the squared absolute value of.
///
/// `ctx` (`anytype`):
/// A context struct providing necessary resources and configuration for the
/// operation. The required fields depend on the output and operand types. If
/// the context is missing required fields or contains unnecessary or wrongly
/// typed fields, the compiler will emit a detailed error message describing the
/// expected structure.
///
/// Returns
/// -------
/// `void`
///
/// Errors
/// ------
/// ``:
///
/// Notes
/// -----
/// When the output and input types are the same, aliasing is allowed.
pub inline fn abs2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isArray(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domain(O)) {
        .expression => @compileError("zml.abs2_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
        .array => switch (comptime types.domain(X)) {
            .array => { // array = abs(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int, .float, .cfloat => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        .integer, .rational, .real, .complex => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        },
                    },
                    .integer, .rational, .real, .complex => switch (types.numericType(types.Numeric(X))) {
                        .bool, .int, .float, .cfloat => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        },
                        .integer, .rational, .real, .complex => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                },
                            );
                        },
                    },
                };

                return array.abs2_(
                    o,
                    x,
                    ctx,
                );
            },
            else => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .matrix => @compileError("zml.abs2_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
        .numeric => switch (comptime types.domain(X)) {
            .numeric => { // numeric = abs(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                int.mul(x, x),
                                .{},
                            ) catch unreachable;
                        },
                        .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                float.mul(x, x),
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                cfloat.abs2(x),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                },
                            );

                            try ops.mul_(
                                o,
                                x,
                                x,
                                ctx,
                            );
                        },
                        .rational => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                },
                            );

                            try ops.mul_(
                                o,
                                x,
                                x,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .integer => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                int.mul(x, x),
                                .{ .allocator = ctx.allocator },
                            );
                        },
                        .float => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                float.mul(x, x),
                                ctx,
                            );
                        },
                        .cfloat => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                cfloat.abs2(x),
                                ctx,
                            );
                        },
                        .integer => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                },
                            );

                            try ops.mul_(
                                o,
                                x,
                                x,
                                ctx,
                            );
                        },
                        .rational => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                },
                            );

                            try ops.mul_(
                                o,
                                x,
                                x,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .rational => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                int.mul(x, x),
                                .{ .allocator = ctx.allocator },
                            );
                        },
                        .float => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                float.mul(x, x),
                                ctx,
                            );
                        },
                        .cfloat => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                cfloat.abs2(x),
                                ctx,
                            );
                        },
                        .integer => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                },
                            );

                            try ops.mul_(
                                o,
                                x,
                                x,
                                ctx,
                            );
                        },
                        .rational => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                },
                            );

                            try ops.mul_(
                                o,
                                x,
                                x,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .real => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .complex => @compileError("zml.abs2_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                }
            },
            else => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.abs2_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
