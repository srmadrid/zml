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

pub fn Abs2(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.abs2: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return cfloat.Abs2(X),
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return complex.Abs2(X),
        .custom => {
            if (comptime !types.hasMethod(X, "Abs2", fn (type) type, &.{X}))
                @compileError("zml.numeric.abs2: " ++ @typeName(X) ++ " must implement `fn Abs2(type) type`");

            return X.Abs2(X);
        },
    }
}

/// Returns the squared absolute value of a numeric `x`. For complex values,
/// returns the squared magnitude.
///
/// ## Signature
/// ```zig
/// numeric.abs2(x: X, ctx: anytype) !numeric.Abs2(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the squared absolute value of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Abs2(X)`.
///
/// #### `numeric.Abs2(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Abs2(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Abs2(@TypeOf(x))`: The squared absolute value of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Abs2(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Abs2` method. The expected signature and
/// behavior of `Abs2` are as follows:
/// * `fn Abs2(type) type`: Returns the return type of `abs2` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Abs2(X)` as `R`. Then, `R` or `X`
/// must implement the required `abs2` method. The expected signatures and
/// behavior of `abs2` are as follows:
/// * `R` is not allocated: `fn abs2(X) R`: Returns the squared absolute value
///   of `x`.
/// * `R` is allocated: `fn abs2(std.mem.Allocator, X) !R`: Returns the
///   squared absolute value of `x` as a newly allocated value.
pub inline fn abs2(x: anytype, ctx: anytype) !numeric.Abs2(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Abs2(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return int.mul(x, x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.mul(x, x);
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return dyadic.mul(x, x);
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.abs2(x);
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

            return integer.mul(ctx.allocator, x, x);
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

            return rational.mul(ctx.allocator, x, x);
        },
        .real => @compileError("zml.numeric.abs2: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.abs2: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "abs2",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.abs2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs2(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.abs2(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "abs2",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.abs2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs2(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.abs2(x);
            }
        },
    }
}
