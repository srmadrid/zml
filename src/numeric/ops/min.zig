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

pub fn Min(X: type, Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.numeric.min: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.haveMethod(
                &.{ X, Y },
                "Min",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.min: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Min(type, type) type`");

            return Impl.Min(X, Y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "Min", fn (type, type) type, &.{ X, Y }))
                @compileError("zml.numeric.min: " ++ @typeName(X) ++ " must implement `fn Min(type, type) type`");

            return X.Min(X, Y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "Min", fn (type, type) type, &.{ X, Y }))
            @compileError("zml.numeric.min: " ++ @typeName(Y) ++ " must implement `fn Min(type, type) type`");

        return Y.Min(X, Y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return bool,
            .int => return int.Min(X, Y),
            .float => return float.Min(X, Y),
            .dyadic => return dyadic.Min(X, Y),
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.Min(X, Y),
            .float => return float.Min(X, Y),
            .dyadic => return dyadic.Min(X, Y),
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Min(X, Y),
            .dyadic => return dyadic.Min(X, Y),
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Min(X, Y),
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int => return integer.Integer,
            .float, .dyadic => return rational.Rational,
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return rational.Rational,
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return real.Real,
            .cfloat => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .integer, .rational, .real => return real.Real,
            .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .custom => unreachable,
        },
        .complex => @compileError("zml.numeric.min: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
        .custom => unreachable,
    }
}

/// Returns the minimum between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.min(x: X, y: Y, ctx: anytype) !numeric.Min(X, Y)
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
/// The fields of `ctx` depend on `numeric.Min(X, Y)`.
///
/// #### `numeric.Min(X, Y)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Min(X, Y)` is allocated and `X != Y`
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// #### `numeric.Min(X, Y)` is allocated and `X == Y`
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// ## Returns
/// `numeric.Min(@TypeOf(x), @TypeOf(y))`: The minimum between `x` and `y`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `min` method. The expected signature
/// and behavior of `Min` are as follows:
/// * `fn Min(type, type) type`: Returns the return type of `min` for the input
///   types.
///
/// Let us denote the return type `numeric.Min(X, Y)` as `R`. Then, `R`, `X` or
/// `Y` must implement the required `min` method. The expected signatures and
/// behavior of `min` are as follows:
/// * `R` is not allocated: `fn min(X, Y) R`: Returns the minimum between `x`
///   and `y`.
/// * `R` is allocated and `X != Y`: `fn min(std.mem.Allocator, X, Y) !R`:
///   Returns the minimum between `x` and `y` as a newly allocated value.
/// * `R` is allocated and `X == Y`: `fn min(?std.mem.Allocator, X, Y) !R`:
///   Returns the minimum between `x` and `y` as a newly allocated value, if the
///   allocator is provided, or as a read-only view if not. If not provided, it
///   must not fail.
pub inline fn min(x: anytype, y: anytype, ctx: anytype) !numeric.Min(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Min(X, Y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            if (comptime types.isAllocated(R)) {
                if (comptime X == Y) {
                    const Impl: type = comptime types.haveMethod(
                        &.{ R, X },
                        "min",
                        fn (?std.mem.Allocator, X, Y) anyerror!R,
                        &.{ ?std.mem.Allocator, X, Y },
                    ) orelse
                        @compileError("zml.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn min(?std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{
                                .type = std.mem.Allocator,
                                .required = false,
                                .description = "The allocator to use for the custom numeric's memory allocation. If not provided, a read-only view will be returned.",
                            },
                        },
                    );

                    return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                        Impl.min(ctx.allocator, x, y)
                    else
                        Impl.min(null, x, y) catch unreachable;
                } else {
                    const Impl: type = comptime types.haveMethod(
                        &.{ R, X, Y },
                        "min",
                        fn (std.mem.Allocator, X, Y) anyerror!R,
                        &.{ std.mem.Allocator, X, Y },
                    ) orelse
                        @compileError("zml.numeric.min: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn min(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                    return Impl.min(ctx.allocator, x, y);
                }
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X, Y },
                    "min",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.min: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn min(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.min(x, y);
            }
        } else { // only X custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "min",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn min(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.min(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "min",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn min(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.min(x, y);
            }
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        if (comptime types.isAllocated(R)) {
            const Impl: type = comptime types.haveMethod(
                &.{ R, Y },
                "min",
                fn (std.mem.Allocator, X, Y) anyerror!R,
                &.{ std.mem.Allocator, X, Y },
            ) orelse
                @compileError("zml.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn min(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

            return Impl.min(ctx.allocator, x, y);
        } else {
            const Impl: type = comptime types.haveMethod(
                &.{ R, Y },
                "min",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.min: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn min(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            comptime types.validateContext(@TypeOf(ctx), .{});

            return Impl.min(x, y);
        }
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => return x and y,
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.min(x, y);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.min(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.min(x, y);
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

                return integer.min(ctx.allocator, x, y);
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

                return rational.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.min(x, y);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.min(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.min(x, y);
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

                return integer.min(ctx.allocator, x, y);
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

                return rational.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.min(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.min(x, y);
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

                return rational.min(ctx.allocator, x, y.asRational());
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

                return rational.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.min(x, y);
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

                return rational.min(ctx.allocator, x, y.asRational());
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

                return rational.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
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

                return integer.min(ctx.allocator, x, y);
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

                return rational.min(ctx.allocator, x.asRational(), y);
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
                    integer.min(ctx.allocator, x, y)
                else
                    integer.min(null, x, y) catch unreachable;
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

                return rational.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
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

                return rational.min(ctx.allocator, x, y);
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

                return rational.min(ctx.allocator, x, y);
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
                    rational.min(ctx.allocator, x, y)
                else
                    rational.min(null, x, y) catch unreachable;
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

                return real.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
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

                return real.min(ctx.allocator, x, y);
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
                    real.min(ctx.allocator, x, y)
                else
                    real.min(null, x, y) catch unreachable;
            },
            .complex => unreachable,
            .custom => unreachable,
        },
        .complex => unreachable,
        .custom => unreachable,
    }
}
