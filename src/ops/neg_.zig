const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const real = @import("../real.zig");
const complex = @import("../complex.zig");

const array = @import("../array.zig");

/// Performs in-place computation of the negation of `x`.
///
/// The `neg_` routine computes the negation of its input `x`, and
/// stores the result directly into `o`, automatically validating the provided
/// context. The operation is performed in the input's precision, and the
/// result is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = neg(Numeric)**: scalar negation
/// - **Numeric = neg(Vector)**: element-wise negation.
/// - **Matrix = neg(Matrix)**: element-wise negation.
/// - **Array = neg(Array)**: element-wise negation.
/// - **Expression = neg(Expression)**: symbolic negation.
///
/// Signature
/// ---------
/// ```zig
/// fn neg_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the negation of.
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
pub inline fn neg_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.neg_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isExpression(O) and !types.isExpression(X) and
        !types.isArray(O) and !types.isArray(X) and
        !types.isMatrix(O) and !types.isMatrix(X) and
        !types.isVector(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domainType(O)) {
        .expression => @compileError("zml.neg_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
        .array => switch (comptime types.domainType(X)) {
            .array => { // array = neg(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        else => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .integer => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
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
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .rational => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
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
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .complex => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            if (types.Scalar(O) == types.Scalar(X)) {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    },
                                );
                            } else {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );
                            }
                        },
                        else => {
                            types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );
                        },
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                };

                return array.neg_(
                    o,
                    x,
                    ctx,
                );
            },
            else => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .matrix => @compileError("zml.neg_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => { // numeric = neg(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                -x,
                                .{},
                            ) catch unreachable;
                        },
                        .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                -x,
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                cfloat.neg(x),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                integer.neg(null, x) catch unreachable,
                                .{},
                            ) catch unreachable;
                        },
                        .rational => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                rational.neg(null, x) catch unreachable,
                                .{},
                            ) catch unreachable;
                        },
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                complex.neg(null, x) catch unreachable,
                                .{},
                            ) catch unreachable;
                        },
                    },
                    .integer => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                -x,
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
                                -x,
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
                                cfloat.neg(x),
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
                                o.positive = !o.positive;
                            } else {
                                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                    try ops.set(
                                        o,
                                        integer.neg(null, x) catch unreachable,
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
                                rational.neg(null, x) catch unreachable,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                complex.neg(null, x) catch unreachable,
                                ctx,
                            );
                        },
                    },
                    .rational => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                -x,
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
                                -x,
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
                                cfloat.neg(x),
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
                                integer.neg(null, x) catch unreachable,
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
                                o.num.positive = !o.num.positive;
                            } else {
                                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                    try ops.set(
                                        o,
                                        rational.neg(null, x) catch unreachable,
                                        .{ .allocator = allocator },
                                    );
                                } else {
                                    return error.AllocatorRequired;
                                }
                            }
                        },
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                complex.neg(null, x) catch unreachable,
                                ctx,
                            );
                        },
                    },
                    .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .complex => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                -x,
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
                                -x,
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
                                cfloat.neg(x),
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
                                integer.neg(null, x) catch unreachable,
                                ctx,
                            );
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
                                rational.neg(null, x) catch unreachable,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.neg_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            if (comptime types.Scalar(O) == types.Scalar(X)) {
                                const spec =
                                    .{
                                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (((comptime types.Scalar(O) == rational.Rational) and o.re.num.limbs == x.re.num.limbs and
                                    o.re.den.limbs == x.re.den.limbs and
                                    o.im.num.limbs == x.im.num.limbs and
                                    o.im.den.limbs == x.im.den.limbs) or
                                    ((comptime types.Scalar(O) == real.Real) and o.re.rational.num.limbs == x.re.rational.num.limbs and
                                        o.re.rational.den.limbs == x.re.rational.den.limbs and
                                        o.im.rational.num.limbs == x.im.rational.num.limbs and
                                        o.im.rational.den.limbs == x.im.rational.den.limbs))
                                {
                                    if (comptime types.Scalar(O) == rational.Rational) {
                                        o.re.num.positive = !o.re.num.positive;
                                        o.im.num.positive = !o.im.num.positive;
                                    } else {
                                        o.re.rational.num.positive = !o.re.rational.num.positive;
                                        o.im.rational.num.positive = !o.im.rational.num.positive;
                                    }
                                } else {
                                    if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                        try ops.set(
                                            o,
                                            complex.neg(null, x) catch unreachable,
                                            .{ .allocator = allocator },
                                        );
                                    } else {
                                        return error.AllocatorRequired;
                                    }
                                }
                            } else {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    complex.neg(null, x) catch unreachable,
                                    ctx,
                                );
                            }
                        },
                    },
                }
            },
            else => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.neg_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
