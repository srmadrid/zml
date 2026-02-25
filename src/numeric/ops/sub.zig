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

pub fn Sub(X: type, Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.numeric.sub: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.haveMethod(
                &.{ X, Y },
                "Sub",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.sub: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Sub(type, type) type`");

            return Impl.Sub(X, Y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "Sub", fn (type, type) type, &.{ X, Y }))
                @compileError("zml.numeric.sub: " ++ @typeName(X) ++ " must implement `fn Sub(type, type) type`");

            return X.Sub(X, Y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "Sub", fn (type, type) type, &.{ X, Y }))
            @compileError("zml.numeric.sub: " ++ @typeName(Y) ++ " must implement `fn Sub(type, type) type`");

        return Y.Sub(X, Y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => @compileError("zml.numeric.sub: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int => return int.Sub(X, Y),
            .float => return float.Sub(X, Y),
            .dyadic => return dyadic.Sub(X, Y),
            .cfloat => return cfloat.Sub(X, Y),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => return int.Sub(X, Y),
            .float => return float.Sub(X, Y),
            .dyadic => return dyadic.Sub(X, Y),
            .cfloat => return cfloat.Sub(X, Y),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Sub(X, Y),
            .dyadic => return dyadic.Sub(X, Y),
            .cfloat => return cfloat.Sub(X, Y),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return dyadic.Sub(X, Y),
            .cfloat => return cfloat.Sub(X, Y),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat => return cfloat.Sub(X, Y),
            .integer, .rational, .real => return complex.Sub(complex.Complex(rational.Rational), Y),
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .integer => switch (comptime types.numericType(Y)) {
            .bool, .int => return integer.Integer,
            .float, .dyadic => return rational.Rational,
            .cfloat => return complex.Sub(complex.Complex(rational.Rational), Y),
            .integer => return integer.Integer,
            .rational => return rational.Rational,
            .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .rational => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return rational.Rational,
            .cfloat => return complex.Sub(complex.Complex(rational.Rational), Y),
            .integer, .rational => return rational.Rational,
            .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .real => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => return real.Real,
            .cfloat => return complex.Sub(complex.Complex(rational.Rational), Y),
            .integer, .rational, .real => return real.Real,
            .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .complex => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real, .complex => return complex.Sub(X, Y),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Performs subtraction between any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.sub(x: X, y: Y, ctx: anytype) !numeric.Sub(X, Y)
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
/// The fields of `ctx` depend on `numeric.Sub(X, Y)`.
///
/// #### `numeric.Sub(X, Y)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Sub(X, Y)` is allocated
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Sub` method. The expected signature
/// and behavior of `Sub` are as follows:
/// * `fn Sub(type, type) type`: Returns the return type of `sub` for the input
///   types.
///
/// Let us denote the return type `numeric.Sub(X, Y)` as `R`. Then, `R`, `X` or
/// `Y` must implement the required `sub` method. The expected signatures and
/// behavior of `sub` are as follows:
/// * `R` is not allocated: `fn sub(X, Y) R`: Returns the subtraction of `x` and
///   `y`.
/// * `R` is allocated: `fn sub(std.mem.Allocator, X, Y) !R`: Returns the
///   sub of `x` and `y` as a newly allocated value.
pub inline fn sub(x: anytype, y: anytype, ctx: anytype) !numeric.Sub(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Sub(X, Y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X, Y },
                    "sub",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.sub: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn sub(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.sub(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X, Y },
                    "sub",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.sub: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn sub(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.sub(x, y);
            }
        } else { // only X custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "sub",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.sub: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sub(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.sub(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "sub",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.sub: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sub(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.sub(x, y);
            }
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        if (comptime types.isAllocated(R)) {
            const Impl: type = comptime types.haveMethod(
                &.{ R, Y },
                "sub",
                fn (std.mem.Allocator, X, Y) anyerror!R,
                &.{ std.mem.Allocator, X, Y },
            ) orelse
                @compileError("zml.numeric.sub: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn sub(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

            return Impl.sub(ctx.allocator, x, y);
        } else {
            const Impl: type = comptime types.haveMethod(
                &.{ R, Y },
                "sub",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.sub: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn sub(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            comptime types.validateContext(@TypeOf(ctx), .{});

            return Impl.sub(x, y);
        }
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => unreachable,
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.sub(x, y);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sub(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.sub(x, y);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sub(x, y);
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

                return integer.sub(ctx.allocator, x, y);
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

                return rational.sub(ctx.allocator, x, y);
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return int.sub(x, y);
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sub(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.sub(x, y);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sub(x, y);
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

                return integer.sub(ctx.allocator, x, y);
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

                return rational.sub(ctx.allocator, x, y);
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.sub(x, y);
            },
            .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.sub(x, y);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sub(x, y);
            },
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

                return rational.sub(ctx.allocator, x, y.asRational());
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

                return rational.sub(ctx.allocator, x, y);
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
        .dyadic => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return dyadic.sub(x, y);
            },
            .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sub(x, y);
            },
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

                return rational.sub(ctx.allocator, x, y.asRational());
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

                return rational.sub(ctx.allocator, x, y);
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
        .cfloat => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return cfloat.sub(x, y);
            },
            .integer, .rational, .real => {
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

                return complex.sub(ctx.allocator, x, y.asComplex());
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
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

                return integer.sub(ctx.allocator, x, y);
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

                return rational.sub(ctx.allocator, x.asRational(), y);
            },
            .cfloat => {
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

                return complex.sub(ctx.allocator, x.asComplex(), y);
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

                return integer.sub(ctx.allocator, x, y);
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

                return rational.sub(ctx.allocator, x, y);
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
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

                return rational.sub(ctx.allocator, x, y);
            },
            .cfloat => {
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

                return complex.sub(ctx.allocator, x.asComplex(), y);
            },
            .integer, .rational => {
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

                return rational.sub(ctx.allocator, x, y);
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
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

                return real.sub(ctx.allocator, x, y);
            },
            .cfloat => {
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

                return complex.sub(ctx.allocator, x.asComplex(), y);
            },
            .integer, .rational, .real => {
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

                return real.sub(ctx.allocator, x, y);
            },
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
        .complex => switch (comptime types.numericType(Y)) {
            .bool, .int, .float, .dyadic, .cfloat, .integer, .rational, .real, .complex => {
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

                return complex.sub(ctx.allocator, x, y);
            },
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}
