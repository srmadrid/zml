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

/// Performs in-place computation of the complex conjugate of `x`.
///
/// The `conj_` routine computes the complex conjugate of its input `x`, and
/// stores the result directly into `o`, automatically validating the provided
/// context. The operation is performed in the input's precision, and the
/// result is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = conj(Numeric)**: scalar complex conjugate.
/// - **Numeric = conj(Vector)**: element-wise complex conjugate.
/// - **Matrix = conj(Matrix)**: element-wise complex conjugate.
/// - **Array = conj(Array)**: element-wise complex conjugate.
/// - **Expression = conj(Expression)**: symbolic complex conjugate.
///
/// Signature
/// ---------
/// ```zig
/// fn conj_(o: *O, x: X, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The operand to compute the complex conjugate of.
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
pub inline fn conj_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.conj_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isExpression(O) and !types.isExpression(X) and
        !types.isArray(O) and !types.isArray(X) and
        !types.isMatrix(O) and !types.isMatrix(X) and
        !types.isVector(X) and
        !types.isNumeric(O) and !types.isNumeric(X))
        @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X));

    switch (comptime types.domainType(O)) {
        .expression => @compileError("zml.conj_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
        .array => switch (comptime types.domainType(X)) {
            .array => { // array = conj(array)
                comptime switch (types.numericType(types.Numeric(O))) {
                    .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        else => {
                            types.validateContext(@TypeOf(ctx), .{});
                        },
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .integer => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
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
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .rational => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
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
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                    .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .complex => switch (types.numericType(types.Numeric(X))) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
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
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    },
                };

                return array.conj_(
                    o,
                    x,
                    ctx,
                );
            },
            else => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        .matrix => @compileError("zml.conj_ between " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " not implemented yet"),
        .numeric => switch (comptime types.domainType(X)) {
            .numeric => { // numeric = conj(numeric)
                switch (comptime types.numericType(O)) {
                    .bool, .int, .float, .cfloat => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                x,
                                .{},
                            ) catch unreachable;
                        },
                        .float => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                x,
                                .{},
                            ) catch unreachable;
                        },
                        .cfloat => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                x.conj(),
                                .{},
                            ) catch unreachable;
                        },
                        .integer => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                x,
                                .{},
                            ) catch unreachable;
                        },
                        .rational => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                x,
                                .{},
                            ) catch unreachable;
                        },
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            comptime types.validateContext(@TypeOf(ctx), .{});

                            ops.set(
                                o,
                                complex.conj(null, x) catch unreachable,
                                .{},
                            ) catch unreachable;
                        },
                    },
                    .integer => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                x,
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
                                x,
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
                                x.conj(),
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
                                // No-op, integer is already its own conjugate
                            } else {
                                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                    try ops.set(
                                        o,
                                        x,
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
                                x,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                complex.conj(null, x) catch unreachable,
                                ctx,
                            );
                        },
                    },
                    .rational => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                x,
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
                                x,
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
                                x.conj(),
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
                                x,
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
                                // No-op, rational is already its own conjugate
                            } else {
                                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                    try ops.set(
                                        o,
                                        x,
                                        .{ .allocator = allocator },
                                    );
                                } else {
                                    return error.AllocatorRequired;
                                }
                            }
                        },
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .complex => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                x,
                                ctx,
                            );
                        },
                    },
                    .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                    .complex => switch (comptime types.numericType(X)) {
                        .bool => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
                        .int => {
                            comptime types.validateContext(
                                @TypeOf(ctx),
                                .{
                                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                                },
                            );

                            try ops.set(
                                o,
                                x,
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
                                x,
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
                                x.conj(),
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
                                x,
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
                                x,
                                ctx,
                            );
                        },
                        .real => @compileError("zml.conj_ not implemented yet for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
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
                                        o.im.num.positive = !o.im.num.positive;
                                    } else {
                                        o.im.rational.num.positive = !o.im.rational.num.positive;
                                    }
                                } else {
                                    if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                                        try ops.set(
                                            o,
                                            complex.conj(null, x) catch unreachable,
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
                                    complex.conj(null, x) catch unreachable,
                                    ctx,
                                );
                            }
                        },
                    },
                }
            },
            else => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
        },
        else => @compileError("zml.conj_ not defined for " ++ @typeName(O) ++ " and " ++ @typeName(X)),
    }
}
