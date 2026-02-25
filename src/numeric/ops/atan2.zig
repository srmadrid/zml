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

pub fn Atan2(X: type, Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y))
        @compileError("zml.numeric.atan2: x and y must be numerics, got \n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            const Impl: type = comptime types.anyHasMethod(
                &.{ X, Y },
                "Atan2",
                fn (type, type) type,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.atan2: " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn Atan2(type, type) type`");

            return Impl.Atan2(X, Y);
        } else { // only X custom
            comptime if (!types.hasMethod(X, "Atan2", fn (type, type) type, &.{ X, Y }))
                @compileError("zml.numeric.atan2: " ++ @typeName(X) ++ " must implement `fn Atan2(type, type) type`");

            return X.Atan2(X, Y);
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        comptime if (!types.hasMethod(Y, "Atan2", fn (type, type) type, &.{ X, Y }))
            @compileError("zml.numeric.atan2: " ++ @typeName(Y) ++ " must implement `fn Atan2(type, type) type`");

        return Y.Atan2(X, Y);
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => @compileError("zml.numeric.atan2: not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ "."),
            .int, .float => return float.Atan2(X, Y),
            .dyadic => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .cfloat => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .integer => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .rational => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .real => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .complex => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .custom => unreachable,
        },
        .int => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Atan2(X, Y),
            .dyadic => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .cfloat => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .integer => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .rational => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .real => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .complex => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .custom => unreachable,
        },
        .float => switch (comptime types.numericType(Y)) {
            .bool, .int, .float => return float.Atan2(X, Y),
            .dyadic => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .cfloat => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .integer, .rational => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .real => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .complex => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
            .custom => unreachable,
        },
        .dyadic => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .cfloat => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .integer => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .rational => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .real => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .complex => @compileError("zml.numeric.atan2: not implemented for " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " yet."),
        .custom => unreachable,
    }
}

/// Computes the arctangent `tan⁻¹(y/x)` of any two numeric operands.
///
/// ## Signature
/// ```zig
/// numeric.atan2(x: X, y: Y, ctx: anytype) !numeric.Atan2(X, Y)
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
/// The fields of `ctx` depend on `numeric.Atan2(X, Y)`.
///
/// #### `numeric.Atan2(X, Y)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Atan2(X, Y)` is allocated
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Atan2(@TypeOf(x), @TypeOf(y))`: The arctangent of `x/y`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Atan2(X, Y)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` or `Y` must implement the required `Atan2` method. The expected
/// signature and behavior of `Atan2` are as follows:
/// * `fn Atan2(type, type) type`: Returns the return type of `atan2` for the
///   input types.
///
/// Let us denote the return type `numeric.Atan2(X, Y)` as `R`. Then, `R`, `X`
/// or `Y` must implement the required `atan2` method. The expected signatures
/// and behavior of `atan2` are as follows:
/// * `R` is not allocated: `fn atan2(X, Y) R`: Returns the arctangent of `x/y`.
/// * `R` is allocated: `fn atan2(std.mem.Allocator, X, Y) !R`: Returns the
///   arctangent of `x/y` as a newly allocated value.
pub inline fn atan2(x: anytype, y: anytype, ctx: anytype) !numeric.Atan2(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = numeric.Atan2(X, Y);

    if (comptime types.isCustomType(X)) {
        if (comptime types.isCustomType(Y)) { // X and Y both custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X, Y },
                    "atan2",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.atan2: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn atan2(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.atan2(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X, Y },
                    "atan2",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.atan2: " ++ @typeName(R) ++ ", " ++ @typeName(X) ++ " or " ++ @typeName(Y) ++ " must implement `fn atan2(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.atan2(x, y);
            }
        } else { // only X custom
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "atan2",
                    fn (std.mem.Allocator, X, Y) anyerror!R,
                    &.{ std.mem.Allocator, X, Y },
                ) orelse
                    @compileError("zml.numeric.atan2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atan2(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.atan2(ctx.allocator, x, y);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "atan2",
                    fn (X, Y) R,
                    &.{ X, Y },
                ) orelse
                    @compileError("zml.numeric.atan2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atan2(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.atan2(x, y);
            }
        }
    } else if (comptime types.isCustomType(Y)) { // only Y custom
        if (comptime types.isAllocated(R)) {
            const Impl: type = comptime types.anyHasMethod(
                &.{ R, Y },
                "atan2",
                fn (std.mem.Allocator, X, Y) anyerror!R,
                &.{ std.mem.Allocator, X, Y },
            ) orelse
                @compileError("zml.numeric.atan2: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn atan2(std.mem.Allocator, " ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") !" ++ @typeName(R) ++ "`");

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

            return Impl.atan2(ctx.allocator, x, y);
        } else {
            const Impl: type = comptime types.anyHasMethod(
                &.{ R, Y },
                "atan2",
                fn (X, Y) R,
                &.{ X, Y },
            ) orelse
                @compileError("zml.numeric.atan2: " ++ @typeName(R) ++ " or " ++ @typeName(Y) ++ " must implement `fn atan2(" ++ @typeName(X) ++ ", " ++ @typeName(Y) ++ ") " ++ @typeName(R) ++ "`");

            comptime types.validateContext(@TypeOf(ctx), .{});

            return Impl.atan2(x, y);
        }
    }

    switch (comptime types.numericType(X)) {
        .bool => switch (comptime types.numericType(Y)) {
            .bool => unreachable,
            .int, .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.atan2(x, y);
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

                return float.atan2(x, y);
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

                return float.atan2(x, y);
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
