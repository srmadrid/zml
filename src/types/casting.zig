const std = @import("std");

const types = @import("../types.zig");
const constants = @import("../constants.zig");
const ops = @import("../ops.zig");

const dyadic = @import("../dyadic.zig");
const integer = @import("../integer.zig");
const Integer = integer.Integer;
const rational = @import("../rational.zig");
const Rational = rational.Rational;
const real = @import("../real.zig");
const Real = real.Real;
const complex = @import("../complex.zig");
const Complex = complex.Complex;

/// Casts a value of any numeric type to any numeric type, without allocation.
/// Some casts may lead to runtime panics if the value cannot be represented
/// in the target type.
///
/// ## Signature
/// ```zig
/// scast(comptime T: type, value: V) T
/// ```
///
/// ## Arguments
/// * `T` (`comptime type`): The type to cast to. Must be a fixed precision
///   numeric type.
/// * `value` (`anytype`): The value to cast.
///
/// ## Returns
/// `T`: The value casted to the type `T`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// ### `T` is a custom numeric type
/// `T` must implement the required `fromV` method, where `V` is the type name,
/// or `Custom` if `V` is also a custom numeric type. The expected signature and
/// behavior of `fromV` are as follows:
/// * `fn fromV(V) T`: Creates a value of type `T` from a value of type `V`.
///
/// ### `V` is a custom numeric type
/// `V` must implement the required `toT` method, where `T` is the type name, or
/// `Custom` if `T` is also a custom numeric type. The expected signature and
/// behavior of `toT` are as follows:
/// * Non-specific types: `fn toT(V, type) T`: Converts a value of type `V` to a
///   value of type `T`.
/// * Specific types: `fn toT(V) T`: Converts a value of type `V` to a value of
///   type `T`.
pub inline fn scast(comptime T: type, value: anytype) T {
    const I: type = @TypeOf(value);
    const O: type = T;

    comptime if (!types.isNumeric(I) or !types.isNumeric(O))
        @compileError("zml.scast: value and T must be numerics, got\n\tvalue: " ++
            @typeName(I) ++ "\n\tT: " ++ @typeName(O) ++ "\n");

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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromBool", fn (bool) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromBool(bool) " ++ @typeName(O) ++ "`");

                return .fromBool(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromInt", fn (I) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromInt(" ++ @typeName(I) ++ ") " ++ @typeName(O) ++ "`");

                return .fromInt(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromFloat", fn (I) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromFloat(" ++ @typeName(I) ++ ") " ++ @typeName(O) ++ "`");

                return .fromFloat(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromDyadic", fn (I) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromDyadic(" ++ @typeName(I) ++ ") " ++ @typeName(O) ++ "`");

                return .fromDyadic(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromCfloat", fn (I) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromCfloat(" ++ @typeName(I) ++ ") " ++ @typeName(O) ++ "`");

                return .fromCfloat(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromInteger", fn (Integer) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromInteger(Integer) " ++ @typeName(O) ++ "`");

                return .fromInteger(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromRational", fn (Rational) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromRational(anytype) " ++ @typeName(O) ++ "`");

                return .fromRational(value);
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
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(O, "fromComplex", fn (I) O, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(O) ++ " must implement `fn fromComplex(" ++ @typeName(I) ++ ") " ++ @typeName(O) ++ "`");

                return .fromComplex(value);
            },
            else => unreachable,
        },
        .custom => switch (comptime types.numericType(O)) {
            .bool => {
                comptime if (!types.hasMethod(I, "toBool", fn (I) bool, &.{@TypeOf(value)}))
                    @compileError("zml.scast: " ++ @typeName(I) ++ " must implement `fn toBool(" ++ @typeName(I) ++ ") bool`");

                return value.toBool();
            },
            .int => {
                comptime if (!types.hasMethod(I, "toInt", fn (I, type) O, &.{ @TypeOf(value), O }))
                    @compileError("zml.scast: " ++ @typeName(I) ++ " must implement `fn toInt(" ++ @typeName(I) ++ ", type) " ++ @typeName(O) ++ "`");

                return value.toInt(O);
            },
            .float => {
                comptime if (!types.hasMethod(I, "toFloat", fn (I, type) O, &.{ @TypeOf(value), O }))
                    @compileError("zml.scast: " ++ @typeName(I) ++ " must implement `fn toFloat(" ++ @typeName(I) ++ ", type) " ++ @typeName(O) ++ "`");

                return value.toFloat(O);
            },
            .dyadic => {
                comptime if (!types.hasMethod(I, "toDyadic", fn (I, type) O, &.{ @TypeOf(value), O }))
                    @compileError("zml.scast: " ++ @typeName(I) ++ " must implement `fn toDyadic(" ++ @typeName(I) ++ ", type) " ++ @typeName(O) ++ "`");

                return value.toDyadic(O);
            },
            .cfloat => {
                comptime if (!types.hasMethod(I, "toCfloat", fn (I, type) O, &.{ @TypeOf(value), O }))
                    @compileError("zml.scast: " ++ @typeName(I) ++ " must implement `fn toCfloat(" ++ @typeName(I) ++ ", type) " ++ @typeName(O) ++ "`");

                return value.toCfloat(O);
            },
            .custom => {
                comptime if (types.isAllocated(O))
                    @compileError("zml.scast: cannot cast to allocated custom type, got\n\tT: " ++ @typeName(O) ++ "\n");

                comptime if (!types.hasMethod(I, "toCustom", fn (I, type) O, &.{ @TypeOf(value), O }))
                    @compileError("zml.scast: " ++ @typeName(I) ++ " must implement `fn toCustom(" ++ @typeName(I) ++ ", type) " ++ @typeName(O) ++ "`");

                return value.toCustom(O);
            },
            else => unreachable,
        },
    }
}

/// Casts a value of any numeric type to any other numeric type. Some casts may
/// lead to runtime panics if the value cannot be represented in the target
/// type.
///
/// ## Signature
/// ```zig
/// cast(comptime T: type, value: V, ctx: anytype) !T
/// ```
///
/// ## Arguments
/// * `T` (`comptime type`): The type to cast to. Must be a supported numeric
///   type.
/// * `value` (`anytype`): The value to cast.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `T` and
///   `V`. If the context is missing required fields or contains unnecessary or
///   wrongly typed fields, the compiler will emit a detailed error message
///   describing the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `T` and `V`.
///
/// #### `T == V` and is allocated
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a view will be returned.
///
/// #### `T` is not allocated
/// The context must be empty.
///
/// #### `T` is allocated
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///   May be optional if a value of type `T` can be created as a view of a value
///   of type `V`.
///
/// ## Returns
/// `T`: The value casted to the type `T`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// ### `T` is a custom numeric type
/// `T` must implement the required `fromV` method, where `V` is the type name,
/// or `Custom` if `V` is also a custom numeric type. The expected signature and
/// behavior of `fromV` are as follows:
/// * Non-allocated: `fn fromV(V) T`: Creates a value of type `T` from a value
///   of type `V`.
/// * Allocated: `fn fromV(std.mem.Allocator, V) !T`: Creates a value of type
///   `T` from a value of type `V` using the provided allocator.
///
/// ## `V` is a custom numeric type
/// `V` must implement the required `toT` method, where `T` is the type name, or
/// `Custom` if `T` is also a custom numeric type. The expected signature and
/// behavior of `toT` are as follows:
/// * Non-specific types:
///   * Non-allocated: `fn toT(V, type) T`: Converts a value of type `V` to a
///     value of type `T`.
///   * Allocated: `fn toT(std.mem.Allocator, V, type) !T`: Converts a value of
///     type `V` to a value of type `T` using the provided allocator.
/// * Specific types:
///   * Non-allocated: `fn toT(V) T`: Converts a value of type `V` to a value of
///     type `T`.
///   * Allocated: `fn toT(std.mem.Allocator, V) !T`: Converts a value of type
///     `V` to a value of type `T` using the provided allocator.
///
/// ## Errors
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value.
pub inline fn cast(
    comptime T: type,
    value: anytype,
    ctx: anytype,
) !T {
    const I: type = @TypeOf(value);
    const O: type = T;

    if (comptime I == O) {
        switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return value;
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the integer's memory allocation. If not provided, a view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    value.copy(ctx.allocator)
                else
                    value.view();
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the rational's memory allocation. If not provided, a view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    value.copy(ctx.allocator)
                else
                    value.view();
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the real's memory allocation. If not provided, a view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    value.copy(ctx.allocator)
                else
                    value.view();
            },
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the complex's memory allocation. If not provided, a view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    value.copy(ctx.allocator)
                else
                    value.view();
            },
            .custom => {
                if (comptime types.isAllocated(O)) {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{
                                .type = std.mem.Allocator,
                                .required = false,
                                .description = "The allocator to use for the custom numeric's memory allocation. If not provided, a view will be returned.",
                            },
                        },
                    );

                    if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator)) {
                        comptime if (!types.hasMethod(O, "copy", fn (O) anyerror!O, &.{O}))
                            @compileError("zml.cast: " ++ @typeName(O) ++ " must implement `fn copy(" ++ @typeName(O) ++ ") !" ++ @typeName(O) ++ "`");

                        return value.copy(ctx.allocator);
                    } else {
                        comptime if (!types.hasMethod(O, "view", fn (O) O, &.{O}))
                            @compileError("zml.cast: " ++ @typeName(O) ++ " must implement `fn view(" ++ @typeName(O) ++ ") " ++ @typeName(O) ++ "`");

                        return value.view();
                    }
                } else {
                    return value;
                }
            },
        }
    }

    switch (comptime types.numericType(I)) {
        .bool => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the integer's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                        },
                    },
                );

                if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator)) {
                    return if (value)
                        constants.one(Integer, ctx)
                    else
                        constants.zero(Integer, ctx);
                } else {
                    return if (value)
                        constants.one(Integer, .{}) catch unreachable
                    else
                        constants.zero(Integer, .{}) catch unreachable;
                }
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the rational's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                        },
                    },
                );

                if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator)) {
                    return if (value)
                        constants.one(Rational, ctx)
                    else
                        constants.zero(Rational, ctx);
                } else {
                    return if (value)
                        constants.one(Rational, .{}) catch unreachable
                    else
                        constants.zero(Rational, .{}) catch unreachable;
                }
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the real's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                        },
                    },
                );

                if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator)) {
                    return if (value)
                        constants.one(Real, ctx)
                    else
                        constants.zero(Real, ctx);
                } else {
                    return if (value)
                        constants.one(Real, .{}) catch unreachable
                    else
                        constants.zero(Real, .{}) catch unreachable;
                }
            },
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the complex's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                        },
                    },
                );

                if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator)) {
                    return if (value)
                        constants.one(O, ctx)
                    else
                        constants.zero(O, ctx);
                } else {
                    return if (value)
                        constants.one(O, .{}) catch unreachable
                    else
                        constants.zero(O, .{}) catch unreachable;
                }
            },
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting bool to allocated custom types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .int => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
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
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                var tv = @import("../int/asRational.zig").asRational(value);
                tv[0].num.limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("zml.cast: casting int to real not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the complex's memory allocation.",
                        },
                    },
                );

                var tv = @import("../int/asComplex.zig").asComplex(value);
                tv[0].re.num.limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting int to allocated custom types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .float => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
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
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                var tv = @import("../float/asRational.zig").asRational(value);
                tv[0].num.limbs = &tv[1][0];
                tv[0].den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("zml.cast: casting float to real not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the complex's memory allocation.",
                        },
                    },
                );

                var tv = @import("../float/asComplex.zig").asComplex(value);
                tv[0].re.num.limbs = &tv[1][0];
                tv[0].re.den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting float to allocated types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .dyadic => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
                    },
                );

                var tv = @import("../dyadic/asInteger.zig").asInteger(value);
                tv[0].limbs = &tv[1];
                return tv[0].copy(ctx.allocator);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                var tv = @import("../dyadic/asRational.zig").asRational(value);
                tv[0].num.limbs = &tv[1][0];
                tv[0].den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("zml.cast: casting dyadic to real not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the complex's memory allocation.",
                        },
                    },
                );

                var tv = @import("../dyadic/asComplex.zig").asComplex(value);
                tv[0].re.num.limbs = &tv[1][0];
                tv[0].re.den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting dyadic to allocated types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .cfloat => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the integer's memory allocation.",
                        },
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
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the rational's memory allocation.",
                        },
                    },
                );

                var tv = @import("../float/asRational.zig").asRational(value.re);
                tv[0].num.limbs = &tv[1][0];
                tv[0].den.limbs = &tv[1][1];
                return tv[0].copy(ctx.allocator);
            },
            .real => @compileError("zml.cast: casting cfloat to real not implemented yet"),
            .complex => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the complex's memory allocation.",
                        },
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
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => unreachable,
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the rational's memory allocation. If not provided, a view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    value.copyToRational(ctx.allocator)
                else
                    value.asRational();
            },
            .real => @compileError("zml.cast: casting integer to real not implemented yet"),
            .complex => @compileError("zml.cast: casting integer to complex not implemented yet"),
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting integer to allocated custom types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .rational => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => @compileError("zml.cast: casting rational to integer not implemented yet"),
            .rational => unreachable,
            .real => @compileError("zml.cast: casting rational to real not implemented yet"),
            .complex => @compileError("zml.cast: casting rational to complex not implemented yet"),
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting rational to allocated custom types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .real => @compileError("zml.cast: casting real to any type not implemented yet"),
        .complex => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer => @compileError("zml.cast: casting complex to integer not implemented yet"),
            .rational => @compileError("zml.cast: casting complex to rational not implemented yet"),
            .real => @compileError("zml.cast: casting complex to real not implemented yet"),
            .complex => @compileError("zml.cast: casting complex to complex not implemented yet"),
            .custom => {
                if (comptime types.isAllocated(O)) {
                    @compileError("zml.cast: casting complex to allocated custom types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
        .custom => switch (comptime types.numericType(O)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return scast(O, value);
            },
            .integer, .rational, .real, .complex => {
                if (comptime types.isAllocated(I)) {
                    @compileError("zml.cast: casting from allocated custom types not implemented yet");
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return scast(O, value);
                }
            },
        },
    }
}
