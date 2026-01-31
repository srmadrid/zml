const std = @import("std");
const types = @import("types.zig");

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");

const _zero: u32 = 0;
const _one: u32 = 1;
const _two: u32 = 2;

/// Returns the additive identity (zero) for the given numeric type `N`.
///
/// If `N` is a custom numeric type, it should implement the required `zero`
/// method. If the method is not found, it must be provided via the context. The
/// expected signature and behavior of `zero` are as follows:
/// * Non-allocated: `fn zero() N`: Returns the zero value.
/// * Allocated: `fn zero(?std.mem.Allocator) !N`: Returns the zero as a newly
///   allocated value, if the allocator is provided, or a read-only view if not.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the zero value for.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `N`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ## Returns
/// * `N`: The zero value.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `N` is an allocated type and an allocator is provided in
///   the context.
///
/// ## Notes
/// * If `N` is an allocated type but no allocator is provided in `ctx`, this
///   returns a read-only view backed by static storage.
///
/// ## Examples
/// * Non-allocated type:
/// ```zig
/// const zero = zml.zero(f64, .{}) catch unreachable;
/// ```
/// * Allocated type (without allocator):
/// ```zig
/// const zero = zml.zero(zml.Integer, .{}) catch unreachable;
/// ```
/// * Allocated type (with allocator):
/// ```zig
/// var zero = try zml.zero(zml.Integer, .{ .allocator = allocator });
/// defer zero.deinit(allocator);
/// ```
/// * Non-allocated custom type (with `zero` method):
/// ```zig
/// const zero = zml.zero(Custom, .{}) catch unreachable;
/// ```
/// * Non-allocated custom type (without `zero` method):
/// ```zig
/// const zero = zml.zero(Custom, .{ .zero = zeroFn }) catch unreachable;
/// ```
/// * Allocated custom type (without allocator, with `zero` method):
/// ```zig
/// var zero = zml.zero(Custom, .{}) catch unreachable;
/// defer zero.deinit(allocator);
/// ```
/// * Allocated custom type (with allocator, with `zero` method):
/// ```zig
/// var zero = try zml.zero(Custom, .{ .allocator = allocator });
/// defer zero.deinit(allocator);
/// ```
/// * Allocated custom type (without allocator, without `zero` method):
/// ```zig
/// var zero = try zml.zero(Custom, .{ .zero = zeroFn }) catch unreachable;
/// defer zero.deinit(allocator);
/// ```
/// * Allocated custom type (with allocator, without `zero` method):
/// ```zig
/// var zero = try zml.zero(Custom, .{ .allocator = allocator, .zero = zeroFn });
/// defer zero.deinit(allocator);
/// ```
pub inline fn zero(
    comptime N: type,
    ctx: anytype,
) !N {
    if (!comptime types.isNumeric(N))
        @compileError("zml.zero: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime types.numericType(N)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return false;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0.0;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .{
                .mantissa = 0,
                .exponent = int.minVal(N.Exponent),
                .positive = true,
            };
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .{
                .re = zero(types.Scalar(N), .{}) catch unreachable,
                .im = zero(types.Scalar(N), .{}) catch unreachable,
            };
        },
        .integer => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the integer's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                return .init(allocator, 2);
            } else {
                return .{
                    .limbs = @ptrCast(@constCast(&_zero)),
                    .size = 0,
                    ._llen = 0,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .rational => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the rational's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var num: integer.Integer = try zero(integer.Integer, .{ .allocator = allocator });
                errdefer num.deinit(allocator);

                return .{
                    .num = num,
                    .den = try one(integer.Integer, .{ .allocator = allocator }),
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .num = zero(integer.Integer, .{}) catch unreachable,
                    .den = one(integer.Integer, .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .real => @compileError("zml.zero: not implemented for " ++ @typeName(N) ++ " yet"),
        .complex => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the complex's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var re: types.Scalar(N) = try zero(types.Scalar(N), .{ .allocator = allocator });
                errdefer re.deinit(allocator);

                return .{
                    .re = re,
                    .im = try zero(types.Scalar(N), .{ .allocator = allocator }),
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .re = zero(types.Scalar(N), .{}) catch unreachable,
                    .im = zero(types.Scalar(N), .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .custom => {
            if (comptime types.isAllocated(N)) {
                if (comptime types.hasMethod(N, "zero", fn (?std.mem.Allocator) anyerror!N, &.{})) {
                    const spec =
                        .{
                            .allocator = .{
                                .type = ?std.mem.Allocator,
                                .required = false,
                                .default = null,
                                .description = "The allocator to use for the custom numeric type's memory allocation. If not provided, a read-only view will be returned.",
                            },
                        };

                    comptime types.validateContext(@TypeOf(ctx), spec);

                    return N.zero(types.getFieldOrDefault(ctx, spec, "allocator"));
                } else {
                    const spec =
                        .{
                            .allocator = .{
                                .type = ?std.mem.Allocator,
                                .required = false,
                                .default = null,
                                .description = "The allocator to use for the custom numeric type's memory allocation. If not provided, a read-only view will be returned.",
                            },
                            .zero = .{
                                .type = fn (allocator: ?std.mem.Allocator) anyerror!N,
                                .required = true,
                                .description = "A function that returns the zero value for the custom numeric type.",
                            },
                        };

                    comptime types.validateContext(@TypeOf(ctx), spec);

                    return ctx.zero(types.getFieldOrDefault(ctx, spec, "allocator"));
                }
            } else {
                if (comptime types.hasMethod(N, "zero", fn () N, &.{})) {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return N.zero();
                } else {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .zero = .{
                                .type = fn () N,
                                .required = true,
                                .description = "A function that returns the zero value for the custom numeric type.",
                            },
                        },
                    );

                    return ctx.zero();
                }
            }
        },
    }
}

/// Returns the multiplicative identity (one) for the given numeric type `N`.
///
/// If `N` is a custom numeric type, it should implement the required `one`
/// method. If the method is not found, it must be provided via the context. The
/// expected signature and behavior of `one` are as follows:
/// * Non-allocated: `fn one() N`: Returns the one value.
/// * Allocated: `fn one(?std.mem.Allocator) !N`: Returns the one as a newly
///   allocated value, if the allocator is provided, or a read-only view if not.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the one value for.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `N`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ## Returns
/// * `N`: The one value.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `N` is an allocated type and an allocator is provided in
///   the context.
///
/// ## Notes
/// * If `N` is an allocated type but no allocator is provided in `ctx`, this
///   returns a read-only view backed by static storage.
///
/// ## Examples
/// * Non-allocated type:
/// ```zig
/// const one = zml.one(f64, .{}) catch unreachable;
/// ```
/// * Allocated type (without allocator):
/// ```zig
/// const one = zml.one(zml.Integer, .{}) catch unreachable;
/// ```
/// * Allocated type (with allocator):
/// ```zig
/// var one = try zml.one(zml.Integer, .{ .allocator = allocator });
/// defer one.deinit(allocator);
/// ```
/// * Non-allocated custom type (with `one` method):
/// ```zig
/// const one = zml.one(Custom, .{}) catch unreachable;
/// ```
/// * Non-allocated custom type (without `one` method):
/// ```zig
/// const one = zml.one(Custom, .{ .one = oneFn }) catch unreachable;
/// ```
/// * Allocated custom type (without allocator, with `one` method):
/// ```zig
/// var one = try zml.one(Custom, .{}) catch unreachable;
/// defer one.deinit(allocator);
/// ```
/// * Allocated custom type (with allocator, with `one` method):
/// ```zig
/// var one = try zml.one(Custom, .{ .allocator = allocator });
/// defer one.deinit(allocator);
/// ```
/// * Allocated custom type (without allocator, without `one` method):
/// ```zig
/// var one = try zml.one(Custom, .{ .one = oneFn }) catch unreachable;
/// defer one.deinit(allocator);
/// ```
/// * Allocated custom type (with allocator, without `one` method):
/// ```zig
/// var one = try zml.one(Custom, .{ .allocator = allocator, .one = oneFn });
/// defer one.deinit(allocator);
/// ```
pub inline fn one(
    comptime N: type,
    ctx: anytype,
) !N {
    if (!comptime types.isNumeric(N))
        @compileError("zml.one: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime types.numericType(N)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return true;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 1;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 1.0;
        },
        .dyadic => @compileError("zml.one: not implemented for dyadic types yet"),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .{
                .re = one(types.Scalar(N), .{}) catch unreachable,
                .im = zero(types.Scalar(N), .{}) catch unreachable,
            };
        },
        .integer => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the integer's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var result: integer.Integer = try .init(allocator, 2);
                result.limbs[0] = 1;
                result.size = 1;
                return result;
            } else {
                return .{
                    .limbs = @ptrCast(@constCast(&_one)),
                    .size = 0,
                    ._llen = 0,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .rational => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the rational's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var num: integer.Integer = try one(integer.Integer, .{ .allocator = allocator });
                errdefer num.deinit(allocator);

                return .{
                    .num = num,
                    .den = try one(integer.Integer, .{ .allocator = allocator }),
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .num = one(integer.Integer, .{}) catch unreachable,
                    .den = one(integer.Integer, .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .real => @compileError("zml.one: not implemented for " ++ @typeName(N) ++ " yet"),
        .complex => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the complex's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var re: types.Scalar(N) = try one(types.Scalar(N), .{ .allocator = allocator });
                errdefer re.deinit(allocator);

                return .{
                    .re = re,
                    .im = try zero(types.Scalar(N), .{ .allocator = allocator }),
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .re = one(types.Scalar(N), .{}) catch unreachable,
                    .im = zero(types.Scalar(N), .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .custom => {
            if (comptime types.isAllocated(N)) {
                if (comptime types.hasMethod(N, "one", fn (?std.mem.Allocator) anyerror!N, &.{})) {
                    const spec =
                        .{
                            .allocator = .{
                                .type = ?std.mem.Allocator,
                                .required = false,
                                .default = null,
                                .description = "The allocator to use for the custom numeric type's memory allocation. If not provided, a read-only view will be returned.",
                            },
                        };

                    comptime types.validateContext(@TypeOf(ctx), spec);

                    return N.one(types.getFieldOrDefault(ctx, spec, "allocator"));
                } else {
                    const spec =
                        .{
                            .allocator = .{
                                .type = ?std.mem.Allocator,
                                .required = false,
                                .default = null,
                                .description = "The allocator to use for the custom numeric type's memory allocation. If not provided, a read-only view will be returned.",
                            },
                            .one = .{
                                .type = fn (allocator: ?std.mem.Allocator) anyerror!N,
                                .required = true,
                                .description = "A function that returns the one value for the custom numeric type.",
                            },
                        };

                    comptime types.validateContext(@TypeOf(ctx), spec);

                    return ctx.one(types.getFieldOrDefault(ctx, spec, "allocator"));
                }
            } else {
                if (comptime types.hasMethod(N, "one", fn () N, &.{})) {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return N.one();
                } else {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .one = .{
                                .type = fn () N,
                                .required = true,
                                .description = "A function that returns the one value for the custom numeric type.",
                            },
                        },
                    );

                    return ctx.one();
                }
            }
        },
    }
}

/// Returns the numeric constant two for the given numeric type `N`.
///
/// If `N` is a custom numeric type, it should implement the required `two`
/// method. If the method is not found, it must be provided via the context. The
/// expected signature and behavior of `two` are as follows:
/// * Non-allocated: `fn two() N`: Returns the two value.
/// * Allocated: `fn two(?std.mem.Allocator) !N`: Returns the two as a newly
///   allocated value, if the allocator is provided, or a read-only view if not.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the two value for.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `N`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ## Returns
/// * `N`: The two value.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `N` is an allocated type and an allocator is provided in
///   the context.
///
/// ## Notes
/// * If `N` is an allocated type but no allocator is provided in `ctx`, this
///   returns a read-only view backed by static storage.
///
/// ## Examples
/// * Non-allocated type:
/// ```zig
/// const two = zml.two(f64, .{}) catch unreachable;
/// ```
/// * Allocated type (without allocator):
/// ```zig
/// const two = zml.two(zml.Integer, .{}) catch unreachable;
/// ```
/// * Allocated type (with allocator):
/// ```zig
/// var two = try zml.two(zml.Integer, .{ .allocator = allocator });
/// defer two.deinit(allocator);
/// ```
/// * Non-allocated custom type (with `two` method):
/// ```zig
/// const two = zml.two(Custom, .{}) catch unreachable;
/// ```
/// * Non-allocated custom type (without `two` method):
/// ```zig
/// const two = zml.two(Custom, .{ .two = twoFn }) catch unreachable;
/// ```
/// * Allocated custom type (without allocator, with `two` method):
/// ```zig
/// var two = try zml.two(Custom, .{}) catch unreachable;
/// defer two.deinit(allocator);
/// ```
/// * Allocated custom type (with allocator, with `two` method):
/// ```zig
/// var two = try zml.two(Custom, .{ .allocator = allocator });
/// defer two.deinit(allocator);
/// ```
/// * Allocated custom type (without allocator, without `two` method):
/// ```zig
/// var two = try zml.two(Custom, .{ .two = twoFn }) catch unreachable;
/// defer two.deinit(allocator);
/// ```
/// * Allocated custom type (with allocator, without `two` method):
/// ```zig
/// var two = try zml.two(Custom, .{ .allocator = allocator, .two = twoFn });
/// defer two.deinit(allocator);
/// ```
pub inline fn two(
    comptime N: type,
    ctx: anytype,
) !N {
    if (!comptime types.isNumeric(N))
        @compileError("zml.two: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime types.numericType(N)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return true;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 2;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 2.0;
        },
        .dyadic => @compileError("zml.two: not implemented for dyadic types yet"),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .{
                .re = two(types.Scalar(N), .{}) catch unreachable,
                .im = zero(types.Scalar(N), .{}) catch unreachable,
            };
        },
        .integer => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the integer's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var result: integer.Integer = try .init(allocator, 2);
                result.limbs[0] = 2;
                result.size = 1;
                return result;
            } else {
                return .{
                    .limbs = @ptrCast(@constCast(&_two)),
                    .size = 0,
                    ._llen = 0,
                    .positive = true,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .rational => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the rational's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var num: integer.Integer = try two(integer.Integer, .{ .allocator = allocator });
                errdefer num.deinit(allocator);

                return .{
                    .num = num,
                    .den = try two(integer.Integer, .{ .allocator = allocator }),
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .num = two(integer.Integer, .{}) catch unreachable,
                    .den = one(integer.Integer, .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .real => @compileError("zml.two: not implemented for " ++ @typeName(N) ++ " yet"),
        .complex => {
            const spec =
                .{
                    .allocator = .{
                        .type = ?std.mem.Allocator,
                        .required = false,
                        .default = null,
                        .description = "The allocator to use for the complex's memory allocation. If not provided, a read-only view backed by static storage will be returned.",
                    },
                };

            comptime types.validateContext(@TypeOf(ctx), spec);

            if (types.getFieldOrDefault(ctx, spec, "allocator")) |allocator| {
                var re: types.Scalar(N) = try two(types.Scalar(N), .{ .allocator = allocator });
                errdefer re.deinit(allocator);

                return .{
                    .re = re,
                    .im = try zero(types.Scalar(N), .{ .allocator = allocator }),
                    .flags = .{ .owns_data = true, .writable = true },
                };
            } else {
                return .{
                    .re = two(types.Scalar(N), .{}) catch unreachable,
                    .im = zero(types.Scalar(N), .{}) catch unreachable,
                    .flags = .{ .owns_data = false, .writable = false },
                };
            }
        },
        .custom => {
            if (comptime types.isAllocated(N)) {
                if (comptime types.hasMethod(N, "two", fn (?std.mem.Allocator) anyerror!N, &.{})) {
                    const spec =
                        .{
                            .allocator = .{
                                .type = ?std.mem.Allocator,
                                .required = false,
                                .default = null,
                                .description = "The allocator to use for the custom numeric type's memory allocation. If not provided, a read-only view will be returned.",
                            },
                        };

                    comptime types.validateContext(@TypeOf(ctx), spec);

                    return N.two(types.getFieldOrDefault(ctx, spec, "allocator"));
                } else {
                    const spec =
                        .{
                            .allocator = .{
                                .type = ?std.mem.Allocator,
                                .required = false,
                                .default = null,
                                .description = "The allocator to use for the custom numeric type's memory allocation. If not provided, a read-only view will be returned.",
                            },
                            .two = .{
                                .type = fn (allocator: ?std.mem.Allocator) anyerror!N,
                                .required = true,
                                .description = "A function that returns the two value for the custom numeric type.",
                            },
                        };

                    comptime types.validateContext(@TypeOf(ctx), spec);

                    return ctx.two(types.getFieldOrDefault(ctx, spec, "allocator"));
                }
            } else {
                if (comptime types.hasMethod(N, "two", fn () N, &.{})) {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return N.two();
                } else {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .two = .{
                                .type = fn () N,
                                .required = true,
                                .description = "A function that returns the two value for the custom numeric type.",
                            },
                        },
                    );

                    return ctx.two();
                }
            }
        },
    }
}
