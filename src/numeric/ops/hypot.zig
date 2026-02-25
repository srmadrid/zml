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

pub fn Hypot(X: type, Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.numeric.hypot: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.anyHasMethod(
                &.{ X, Y },
                "Hypot",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.hypot: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Hypot(type, type) type`");

            return Impl.Hypot(X, Y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "Hypot", fn (type, type) type, &.{ X, Y }))
                @compileError("zml.numeric.hypot: " ++ @typeName(X) ++ " must implement `fn Hypot(type, type) type`");

            return X.Hypot(X, Y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "Hypot", fn (type, type) type, &.{ X, Y }))
            @compileError("zml.numeric.hypot: " ++ @typeName(Y) ++ " must implement `fn Hypot(type, type) type`");

        return Y.Hypot(X, Y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => @compileError("zml.numeric.hypot: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int, .float => return float.Hypot(X, Y),
            .dyadic => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .cfloat => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .integer => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .rational => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .real => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .complex => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Hypot(X, Y),
            .dyadic => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .cfloat => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .integer => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .rational => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .real => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .complex => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Hypot(X, Y),
            .dyadic => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .cfloat => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .integer, .rational => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .real => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .complex => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .custom => unreachable,
        },
        .dyadic => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .cfloat => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .integer => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .rational => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .real => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .complex => @compileError("zml.numeric.hypot: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .custom => unreachable,
    }
}

/// Computes the hypotenuse `√(x² + y²)` of any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.hypot(x: X, y: Y, ctx: anytype) !numeric.Hypot(X, Y)
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
/// The fields of `ctx` depend on `numeric.Hypot(X, Y)`.
///
/// #### `numeric.Hypot(X, Y)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Hypot(X, Y)` is allocated
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Hypot(@TypeOf(x), @TypeOf(y))`: The hypotenuse of `x` and `y`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Hypot(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Hypot` method. The expected
/// signature and behavior of `Hypot` are as follows:
/// * `fn Hypot(type, type) type`: Returns the return type of `hypot` for the
///   input types.
///
/// Let us denote the return type `numeric.Hypot(X, Y)` as `R`. Then, `R`, `X`
/// or `Y` must implement the required `hypot` method. The expected signatures
/// and behavior of `hypot` are as follows:
/// * `R` is not allocated: `fn hypot(X, Y) R`: Returns the hypotenuse of `x`
///   and `y`.
/// * `R` is allocated: `fn hypot(std.mem.Allocator, X, Y) !R`: Returns the
///   hypotenuse of `x` and `y` as a newly allocated value.
pub inline fn hypot(x: anytype, y: anytype, ctx: anytype) !numeric.Hypot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Hypot(X, Y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X, Y },
                    "hypot",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.hypot: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn hypot(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.hypot(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X, Y },
                    "hypot",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.hypot: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn hypot(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.hypot(x, y);
            }
        } else { // only X custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "hypot",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.hypot: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn hypot(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.hypot(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "hypot",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.hypot: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn hypot(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.hypot(x, y);
            }
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        if (comptime types.isAllocated(R)) {
            const Impl: type = comptime types.anyHasMethod(
                &.{ R, Y },
                "hypot",
                fn (std.mem.Allocator, X, Y) anyerror!R,
                &.{ std.mem.Allocator, X, Y },
            ) orelse
                @compileError("zml.numeric.hypot: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn hypot(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

            return Impl.hypot(ctx.allocator, x, y);
        } else {
            const Impl: type = comptime types.anyHasMethod(
                &.{ R, Y },
                "hypot",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.hypot: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn hypot(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            comptime types.validateContext(@TypeOf(ctx), .{});

            return Impl.hypot(x, y);
        }
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => unreachable,
            .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.hypot(x, y);
            },
            .dyadic => unreachable,
            .cfloat => unreachable,
            .integer => unreachable,
            .rational => unreachable,
            .real => unreachable,
            .complex => unreachable,
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.hypot(x, y);
            },
            .dyadic => unreachable,
            .cfloat => unreachable,
            .integer => unreachable,
            .rational => unreachable,
            .real => unreachable,
            .complex => unreachable,
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.hypot(x, y);
            },
            .dyadic => unreachable,
            .cfloat => unreachable,
            .integer => unreachable,
            .rational => unreachable,
            .real => unreachable,
            .complex => unreachable,
            .custom => unreachable,
        },
        .dyadic => unreachable,
        .cfloat => unreachable,
        .integer => unreachable,
        .rational => unreachable,
        .real => unreachable,
        .complex => unreachable,
        .custom => unreachable,
    }
}
