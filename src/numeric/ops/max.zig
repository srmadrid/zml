const std = @import("std");

const types = @import("../../types.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const cfloat = @import("../../cfloat.zig");
const integer = @import("../../integer.zig");
const rational = @import("../../rational.zig");
const real = @import("../../real.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Max(X: type, Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.numeric.max: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            if (comptime !types.hasMethod(X, "Max", fn (type, type) type, &.{ X, Y })) {
                if (comptime !types.hasMethod(Y, "Max", fn (type, type) type, &.{ X, Y }))
                    @compileError("zml.numeric.max: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Max(type, type) type`");

                return Y.Max(X, Y);
            } else {
                return X.Max(X, Y);
            }
        } else { // only X custom
            comptime if (!types.hasMethod(X, "Max", fn (type, type) type, &.{ X, Y }))
                @compileError("zml.numeric.max: " ++ @typeName(X) ++ " must implement `fn Max(type, type) type`");

            return X.Max(X, Y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "Max", fn (type, type) type, &.{ X, Y }))
            @compileError("zml.numeric.max: " ++ @typeName(Y) ++ " must implement `fn Max(type, type) type`");

        return Y.Max(X, Y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return bool,
            .int => return int.Max(X, Y),
            .float => return float.Max(X, Y),
            .dyadic => return dyadic.Max(X, Y),
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.Max(X, Y),
            .float => return float.Max(X, Y),
            .dyadic => return dyadic.Max(X, Y),
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Max(X, Y),
            .dyadic => return dyadic.Max(X, Y),
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Max(X, Y),
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int => return integer.Integer,
            .float, .dyadic => return rational.Rational,
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return rational.Rational,
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return real.Real,
            .cfloat => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational, .real => return real.Real,
            .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zml.numeric.max: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}

/// Returns the maximum between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.max(x: X, y: Y, ctx: anytype) !numeric.Max(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X` and
///   `Y`. If the context is missing required fields or contains unnecessary or
///   wrongly typed fields, the compiler will emit a detailed error message
///   describing the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Max(X, Y)`.
///
/// #### `numeric.Max(X, Y)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Max(X, Y)` is allocated and `X != Y`
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// #### `numeric.Max(X, Y)` is allocated and `X == Y`
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// ## Returns
/// `numeric.Max(@TypeOf(x), @TypeOf(y))`: The maximum between `x` and `y`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// ### `X` is a custom numeric type
/// `X` must implement the required `Max` and `max` methods. The expected
/// signature and behavior of `Max` and `max` are as follows:
/// * `fn Max(type) type`: Returns the return type of `max` for the custom
///   numeric type.
/// * `X.Max(X, Y)` is non-allocated: `fn max(X, Y) X.Max(X, Y)`: Returns the
///   maximum between `x` and `y`.
/// * `X.Max(X, Y)` is allocated and `X != Y`: `fn max(std.mem.Allocator, X, Y) !numeric.Max(X, Y)`:
///   Returns the maximum between `x` and `y` as a newly allocated value.
/// * `X.Max(X, Y)` is allocated and `X == Y`: `fn max(?std.mem.Allocator, X, Y) !numeric.Max(X, Y)`:
///   Returns the maximum between `x` and `y` as a newly allocated value, if the
///   allocator is provided, or as a read-only view if not.
///
/// ### `Y` is a custom numeric type
/// `Y` must implement the required `Max` and `max` methods, with the same
/// specifications as above.
///
/// ### Both `X` and `Y` are custom numeric types
/// At least one of `X` and `Y` must implement the required `Max` and `max`
/// methods, with the same specifications as above. If both implement them,
/// `X`'s implementation will be used.
pub inline fn max(x: anytype, y: anytype, ctx: anytype) !numeric.Max(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Max(X, Y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            if (comptime types.isAllocated(R)) {
                if (comptime !types.hasMethod(X, "max", fn (std.mem.Allocator, X, Y) anyerror!R, &.{ std.mem.Allocator, X, Y })) {
                    comptime if (!types.hasMethod(Y, "max", fn (std.mem.Allocator, X, Y) anyerror!R, &.{ std.mem.Allocator, X, Y }))
                        @compileError("zml.numeric.max: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn max(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{
                                .type = std.mem.Allocator,
                                .required = true,
                                .description = "The allocator to use for the custom numeric's memory allocation.",
                            },
                        },
                    );

                    return Y.max(ctx.allocator, x, y);
                } else {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{
                                .type = std.mem.Allocator,
                                .required = true,
                                .description = "The allocator to use for the custom numeric's memory allocation.",
                            },
                        },
                    );

                    return X.max(ctx.allocator, x, y);
                }
            } else {
                if (comptime !types.hasMethod(X, "max", fn (X, Y) R, &.{ X, Y })) {
                    comptime if (!types.hasMethod(Y, "max", fn (X, Y) R, &.{ X, Y }))
                        @compileError("zml.numeric.max: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn max(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return Y.max(x, y);
                } else {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    return X.max(x, y);
                }
            }
        } else { // only X custom
            if (comptime types.isAllocated(R)) {
                comptime if (!types.hasMethod(X, "max", fn (std.mem.Allocator, X, Y) anyerror!R, &.{ std.mem.Allocator, X, Y }))
                    @compileError("zml.numeric.max: " ++ @typeName(X) ++ " must implement `fn max(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the custom numeric's memory allocation.",
                        },
                    },
                );

                return X.max(ctx.allocator, x, y);
            } else {
                comptime if (!types.hasMethod(X, "max", fn (X, Y) R, &.{ X, Y }))
                    @compileError("zml.numeric.max: " ++ @typeName(Y) ++ " must implement `fn max(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.max(x, y);
            }
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        if (comptime types.isAllocated(R)) {
            comptime if (!types.hasMethod(Y, "max", fn (std.mem.Allocator, X, Y) anyerror!R, &.{ std.mem.Allocator, X, Y }))
                @compileError("zml.numeric.max: " ++ @typeName(Y) ++ " must implement `fn max(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = true,
                        .description = "The allocator to use for the custom numeric's memory allocation.",
                    },
                },
            );

            return Y.max(ctx.allocator, x, y);
        } else {
            comptime if (!types.hasMethod(Y, "max", fn (X, Y) R, &.{ X, Y }))
                @compileError("zml.numeric.max: " ++ @typeName(Y) ++ " must implement `fn max(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            comptime types.validateContext(@TypeOf(ctx), .{});

            return Y.max(x, y);
        }
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return x or y,
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.max(x, y);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.max(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.max(x, y);
            },
            .cfloat => unreachable,
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

                return integer.max(ctx.allocator, x, y);
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

                return rational.max(ctx.allocator, x, y);
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.max(x, y);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.max(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.max(x, y);
            },
            .cfloat => unreachable,
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

                return integer.max(ctx.allocator, x, y);
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

                return rational.max(ctx.allocator, x, y);
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.max(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.max(x, y);
            },
            .cfloat => unreachable,
            .integer => {
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

                return rational.max(ctx.allocator, x, y.asRational());
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

                return rational.max(ctx.allocator, x, y);
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.max(x, y);
            },
            .cfloat => unreachable,
            .integer => {
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

                return rational.max(ctx.allocator, x, y.asRational());
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

                return rational.max(ctx.allocator, x, y);
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .cfloat => unreachable,
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int => {
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

                return integer.max(ctx.allocator, x, y);
            },
            .float, .dyadic => {
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

                return rational.max(ctx.allocator, x.asRational(), y);
            },
            .cfloat => unreachable,
            .integer => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the integer's memory allocation. If not provided, a read-only view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    integer.max(ctx.allocator, x, y)
                else
                    integer.max(null, x, y) catch unreachable;
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

                return rational.max(ctx.allocator, x, y);
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => {
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

                return rational.max(ctx.allocator, x, y);
            },
            .cfloat => unreachable,
            .integer => {
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

                return rational.max(ctx.allocator, x, y);
            },
            .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the rational's memory allocation. If not provided, a read-only view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    rational.max(ctx.allocator, x, y)
                else
                    rational.max(null, x, y) catch unreachable;
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .cfloat => unreachable,
            .integer, .rational => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = true,
                            .description = "The allocator to use for the real's memory allocation.",
                        },
                    },
                );

                return real.max(ctx.allocator, x, y);
            },
            .real => {
                comptime types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{
                            .type = std.mem.Allocator,
                            .required = false,
                            .description = "The allocator to use for the real's memory allocation. If not provided, a read-only view will be returned.",
                        },
                    },
                );

                return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                    real.max(ctx.allocator, x, y)
                else
                    real.max(null, x, y) catch unreachable;
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .complex => unreachable,
        .custom => unreachable,
    }
}
