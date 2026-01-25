const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");

const dyadic = @import("../dyadic.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;
const rational = @import("../rational.zig");
const Rational = rational.Rational;
const real = @import("../real.zig");
const Real = real.Real;
const complex = @import("../complex.zig");
const Complex = complex.Complex;

/// Casts a value of any numeric type to any fixed precision numeric type.
///
/// It is a more concise version of `cast` that does not require any options and
/// cannot return an error.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to cast to. Must be a fixed precision numeric
/// type.
///
/// value (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`): The value to cast.
///
/// Returns
/// -------
/// `T`: The value casted to the type `T`.
///
/// Notes
/// -----
/// This function does not check if the cast is safe, which could lead to
/// runtime panics if the value cannot be represented in the target type.
pub inline fn scast(
    comptime T: type,
    value: anytype,
) T {
    const I: type = @TypeOf(value);
    const O: type = T;

    comptime if (types.isAllocated(O))
        @compileError("Expected a fixed precision type, but got " ++ @typeName(O));

    if (comptime I == O)
        return value;

    switch (comptime types.numericType(I)) {
        .bool => switch (comptime types.numericType(O)) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1.0 else 0.0,
            .dyadic => return .initSet(value),
            .cfloat => return .{
                .re = if (value) 1.0 else 0.0,
                .im = 0.0,
            },
            else => unreachable,
        },
        .int => switch (comptime types.numericType(O)) {
            .bool => return value != 0,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .dyadic => return .initSet(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0.0,
            },
            else => unreachable,
        },
        .float => switch (comptime types.numericType(O)) {
            .bool => return value != 0.0,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .dyadic => return .initSet(value),
            .cfloat => return if (comptime I == types.Scalar(O)) .{
                .re = value,
                .im = 0.0,
            } else .{
                .re = @floatCast(value),
                .im = 0.0,
            },
            else => unreachable,
        },
        .dyadic => switch (comptime types.numericType(O)) {
            .bool => return dyadic.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .dyadic => return .initSet(value),
            .cfloat => return .{
                .re = value.toFloat(types.Scalar(O)),
                .im = 0.0,
            },
            else => unreachable,
        },
        .cfloat => switch (comptime types.numericType(O)) {
            .bool => return value.re != 0 or value.im != 0,
            .int => return @intFromFloat(value.re),
            .float => return if (comptime types.Scalar(I) == O) value.re else @floatCast(value.re),
            .dyadic => return .initSet(value.re),
            .cfloat => return .{
                .re = @floatCast(value.re),
                .im = @floatCast(value.im),
            },
            else => unreachable,
        },
        .integer => switch (comptime types.numericType(O)) {
            .bool => return integer.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .dyadic => return .initSet(value),
            .cfloat => return .{
                .re = value.toFloat(types.Scalar(O)),
                .im = 0.0,
            },
            else => unreachable,
        },
        .rational => switch (comptime types.numericType(O)) {
            .bool => return rational.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .dyadic => return .initSet(value),
            .cfloat => return .{
                .re = value.toFloat(types.Scalar(O)),
                .im = 0.0,
            },
            else => unreachable,
        },
        .real => @compileError("Not implemented yet."),
        .complex => switch (comptime types.numericType(O)) {
            .bool => return complex.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .dyadic => return .initSet(value.re),
            .cfloat => return value.toCFloat(O),
            else => unreachable,
        },
        .custom => unreachable,
    }
}

/// Casts a value of any numeric type to any other numeric type.
///
/// Parameters
/// ----------
/// comptime `T` (`type`): The type to cast to. Must be a supported numeric type.
///
/// `value` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex`):The value to cast.
///
/// `ctx` (`struct`): The context for the cast operation.
///
/// Returns
/// -------
/// `T`: The value casted to the type `T`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value.
///
/// Notes
/// -----
/// This function does not check if the cast is safe, which could lead to
/// runtime panics.
pub inline fn cast(
    comptime T: type,
    value: anytype,
    ctx: anytype,
) !T {
    const I: type = @TypeOf(value);
    const O: type = T;

    if (I == O) {
        switch (comptime types.numericType(O)) {
            .bool, .int, .float, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value;
            },
            .integer => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
            .rational => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
            .real => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
            .complex => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return value.copy(allocator);
                } else {
                    return value;
                }
            },
        }
        return value;
    }

    switch (comptime types.numericType(I)) {
        .bool => switch (comptime types.numericType(O)) {
            .bool => unreachable,
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return if (value) 1 else 0;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return if (value) 1.0 else 0.0;
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .initSet(value);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .{
                    .re = if (value) 1.0 else 0.0,
                    .im = 0.0,
                };
            },
            .integer => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(Integer, .{ .allocator = allocator })
                    else
                        constants.zero(Integer, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(Integer, .{}) catch unreachable
                    else
                        constants.zero(Integer, .{}) catch unreachable;
                }
            },
            .rational => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(Rational, .{ .allocator = allocator })
                    else
                        constants.zero(Rational, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(Rational, .{}) catch unreachable
                    else
                        constants.zero(Rational, .{}) catch unreachable;
                }
            },
            .real => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(Real, .{ .allocator = allocator })
                    else
                        constants.zero(Real, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(Real, .{}) catch unreachable
                    else
                        constants.zero(Real, .{}) catch unreachable;
                }
            },
            .complex => {
                const spec =
                    .{
                        .allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                    };

                comptime types.validateContext(@TypeOf(ctx), spec);

                if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                    return if (value)
                        constants.one(O, .{ .allocator = allocator })
                    else
                        constants.zero(O, .{ .allocator = allocator });
                } else {
                    return if (value)
                        constants.one(O, .{}) catch unreachable
                    else
                        constants.zero(O, .{}) catch unreachable;
                }
            },
        },
        .int => switch (comptime types.numericType(O)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value != 0;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return @intCast(value);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return @floatFromInt(value);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .initSet(value);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .{
                    .re = @floatFromInt(value),
                    .im = 0.0,
                };
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../int/asInteger.zig").asInteger(value);
                tv[0].limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../int/asRational.zig").asRational(value);
                tv[0].num.limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("Not implemented yet: casting from int to real"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../int/asComplex.zig").asComplex(value);
                tv[0].re.num.limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
        },
        .float => switch (comptime types.numericType(O)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value != 0.0;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return @intFromFloat(value);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return @floatCast(value);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .initSet(value);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return if (comptime I == types.Scalar(O)) .{
                    .re = value,
                    .im = 0.0,
                } else .{
                    .re = @floatCast(value),
                    .im = 0.0,
                };
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../float/asInteger.zig").asInteger(value);
                tv[0].limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../float/asRational.zig").asRational(value);
                tv[0].num.limbs = &tv[1][0];
                tv[0].den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("Not implemented yet: casting from float to real"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../float/asComplex.zig").asComplex(value);
                tv[0].re.num.limbs = &tv[1][0];
                tv[0].re.den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
        },
        .dyadic => switch (comptime types.numericType(O)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.ne(value, 0);
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value.toInt(O);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value.toFloat(O);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .initSet(value);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .{
                    .re = value.toFloat(types.Scalar(O)),
                    .im = 0.0,
                };
            },
            .integer => @compileError("Not implemented yet: casting from dyadic to integer"),
            .rational => @compileError("Not implemented yet: casting from dyadic to rational"),
            .real => @compileError("Not implemented yet: casting from dyadic to real"),
            .complex => @compileError("Not implemented yet: casting from dyadic to complex"),
        },
        .cfloat => switch (comptime types.numericType(O)) {
            .bool => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value.re != 0.0 or value.im != 0.0;
            },
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return @intFromFloat(value.re);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return if (comptime types.Scalar(I) == O)
                    value.re
                else
                    @floatCast(value.re);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return .initSet(value.re);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                if (comptime I == O) {
                    return value;
                } else {
                    return .{
                        .re = @floatCast(value.re),
                        .im = @floatCast(value.im),
                    };
                }
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../float/asInteger.zig").asInteger(value.re);
                tv[0].limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../float/asRational.zig").asRational(value.re);
                tv[0].num.limbs = &tv[1][0];
                tv[0].den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("Not implemented yet: casting from cfloat to real"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                var tv = @import("../cfloat/asComplex.zig").asComplex(value);
                tv[0].re.num.limbs = &tv[1][0];
                tv[0].re.den.limbs = &tv[1][1];
                tv[0].im.num.limbs = &tv[1][2];
                tv[0].im.den.limbs = &tv[1][3];
                return tv[0].copy(ctx.allocator);
            },
        },
        .integer => switch (comptime types.numericType(O)) {
            .bool => return integer.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .cfloat => return .{
                .re = value.toFloat(types.Scalar(O)),
                .im = 0.0,
            },
            .integer => unreachable,
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return value.copyToRational(ctx.allocator);
            },
            .real => @compileError("Not implemented yet: casting from integer to real"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return value.copyToComplex(ctx.allocator);
            },
        },
        .rational => switch (comptime types.numericType(O)) {
            .bool => return rational.ne(value, 0),
            .int => return value.toInt(O),
            .float => return value.toFloat(O),
            .cfloat => return .{
                .re = value.toFloat(types.Scalar(O)),
                .im = 0.0,
            },
            .integer => @compileError("Not implemented yet: casting from rational to integer"),
            .rational => unreachable,
            .real => @compileError("Not implemented yet: casting from rational to real"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );

                return value.copyToComplex(ctx.allocator);
            },
        },
        .real => switch (comptime types.numericType(O)) {
            .bool => @compileError("Not implemented yet: casting from real to bool"),
            .int => @compileError("Not implemented yet: casting from real to int"),
            .float => @compileError("Not implemented yet: casting from real to float"),
            .cfloat => @compileError("Not implemented yet: casting from real to cfloat"),
            .integer => @compileError("Not implemented yet: casting from real to integer"),
            .rational => @compileError("Not implemented yet: casting from real to rational"),
            .real => unreachable,
            .complex => @compileError("Not implemented yet: casting from real to complex"),
        },
        .complex => switch (comptime types.numericType(O)) {
            .bool => return complex.ne(value, 0),
            .int => return value.re.toInt(O),
            .float => return value.re.toFloat(O),
            .cfloat => return .{
                .re = value.re.toFloat(types.Scalar(O)),
                .im = value.im.toFloat(types.Scalar(O)),
            },
            .integer => @compileError("Not implemented yet: casting from complex to integer"),
            .rational => return value.re.copy(ctx.allocator),
            .real => @compileError("Not implemented yet: casting from complex to real"),
            .complex => @compileError("Not implemented yet: casting from complex to complex"),
        },
    }
}
