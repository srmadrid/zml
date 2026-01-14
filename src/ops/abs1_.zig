const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");

const array = @import("../array.zig");

/// Performs in-place computation of the absolute value of `x` for real inputs,
/// or the sum of the absolute values of the real and imaginary parts for
/// complex inputs.
///
/// The `abs1_` routine computes the absolute value of its input `x`, and stores
/// the result directly into `o`, automatically validating the provided
/// context. The operation is performed in the input's precision, and the
/// result is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = abs1(Numeric)**: scalar absolute value.
/// - **Array = abs1(Array)**: element-wise absolute value.
/// - **Expression = abs1(Expression)**: symbolic absolute value.
///
/// Signature
/// ---------
/// ```zig
/// fn abs1_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the absolute value of.
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
///
/// When the input and output types are the same, the context may mark the
/// allocator as optional. In practice, the allocator is only truly optional if
/// aliasing is detected. If there is no aliasing, an allocator is still
/// required even though the context marks it as optional.
pub inline fn abs1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs1_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isExpression(O) and !types.isExpression(X) and
        !types.isArray(O) and !types.isArray(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domain(O)) {
        .array => switch (comptime types.domain(X)) {
            .array, .numeric => { // array = abs1(numeric), array = abs1(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        else => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        .real, .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .integer => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .integer => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                },
                            );
                        },
                        else => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        },
                        .real, .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .rational => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .rational => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                },
                            );
                        },
                        else => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        },
                        .real, .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .real, .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                };

                return array.abs1_(
                    o,
                    x,
                    types.stripStruct(ctx, &.{"array_allocator"}),
                );
            },
            else => unreachable,
        },
        .numeric => switch (comptime types.domain(X)) {
            .numeric => { // numeric = abs1(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                int.abs(x),
                                .{},
                            ) catch unreachable;
                        },
                        .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                float.abs(x),
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                cfloat.abs1(x),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                integer.abs(null, x) catch unreachable,
                                .{},
                            ) catch unreachable;
                        },
                        .rational => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                rational.abs(null, x) catch unreachable,
                                .{},
                            ) catch unreachable;
                        },
                        .real => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .integer => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                int.abs(x),
                                ctx,
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
                                float.abs(x),
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
                                cfloat.abs(x),
                                ctx,
                            );
                        },
                        .integer => {
                            const spec =
                                .{
                                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                };

                            comptime types.validateContext(@TypeOf(ctx), spec);

                            if (o.limbs == x.limbs) {
                                o.positive = true;
                            } else {
                                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                    try ops.set(
                                        o,
                                        integer.abs(null, x) catch unreachable,
                                        .{ .allocator = allocator },
                                    );
                                } else {
                                    return error.AllocatorRequired;
                                }
                            }
                        },
                        .rational => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                rational.abs(null, x) catch unreachable,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .rational => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                int.abs(x),
                                ctx,
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
                                float.abs(x),
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
                                cfloat.abs(x),
                                ctx,
                            );
                        },
                        .integer => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                integer.abs(null, x) catch unreachable,
                                ctx,
                            );
                        },
                        .rational => {
                            const spec =
                                .{
                                    .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                };

                            comptime types.validateContext(@TypeOf(ctx), spec);

                            if (o.num.limbs == x.num.limbs and
                                o.den.limbs == x.den.limbs)
                            {
                                o.num.positive = true;
                            } else {
                                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                    try ops.set(
                                        o,
                                        rational.abs(null, x) catch unreachable,
                                        .{ .allocator = allocator },
                                    );
                                } else {
                                    return error.AllocatorRequired;
                                }
                            }
                        },
                        .real => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .real => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .complex => @compileError("zml.abs1_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                }
            },
            else => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.abs1_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
