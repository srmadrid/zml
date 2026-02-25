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

pub fn Sin(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.sin: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Sin(X),
        .int => return float.Sin(X),
        .float => return float.Sin(X),
        .dyadic => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Sin", fn (type) type, &.{X}))
                @compileError("zml.numeric.sin: " ++ @typeName(X) ++ " must implement `fn Sin(type) type`");

            return X.Sin(X);
        },
    }
}

/// Returns the sine `sin(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.sin(x: X, ctx: anytype) !numeric.Sin(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the sine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Sin(X)`.
///
/// #### `numeric.Sin(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Sin(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Sin(@TypeOf(x))`: The sine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Sin(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sin` method. The expected signature and
/// behavior of `Sin` are as follows:
/// * `fn Sin(type) type`: Returns the return type of `sin` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Sin(X)` as `R`. Then, `R` or `X`
/// must implement the required `sin` method. The expected signatures and
/// behavior of `sin` are as follows:
/// * `R` is not allocated: `fn sin(X) R`: Returns the sine of `x`.
/// * `R` is allocated: `fn sin(std.mem.Allocator, X) !R`: Returns the sine of
///   `x` as a newly allocated value.
pub inline fn sin(x: anytype, ctx: anytype) !numeric.Sin(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Sin(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sin(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sin(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sin(x);
        },
        .dyadic => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.sin(x);
        },
        .integer => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sin: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "sin",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.sin: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sin(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.sin(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "sin",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.sin: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sin(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.sin(x);
            }
        },
    }
}
