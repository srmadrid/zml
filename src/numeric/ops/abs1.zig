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

pub fn Abs1(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.abs1: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return cfloat.Abs1(X),
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return complex.Abs1(X),
        .custom => {
            if (comptime !types.hasMethod(X, "Abs1", fn (type) type, &.{X}))
                @compileError("zml.numeric.abs1: " ++ @typeName(X) ++ " must implement `fn Abs1(type) type`");

            return X.Abs1(X);
        },
    }
}

/// Returns the 1-norm of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.abs1(x: X, ctx: anytype) !numeric.Abs1(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the 1-norm of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Abs1(X)` and `X`.
///
/// #### `numeric.Abs1(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Abs1(X)` is allocated and `X.has_simple_abs1` exists and is true
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// #### `numeric.Abs1(X)` is allocated and `X.has_simple_abs1` does not exist or is false
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Abs1(@TypeOf(x))`: The 1-norm of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Abs1(X)` is allocated and an allocator is
///   provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Abs1` method. The expected signature and
/// behavior of `Abs1` are as follows:
/// * `fn Abs1(type) type`: Returns the return type of `abs1` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Abs1(X)` as `R`. Then, `R` or `X`
/// must implement the required `abs1` method. The expected signatures and
/// behavior of `abs1` are as follows:
/// * `R` is not allocated: `fn abs1(X) R`: Returns the 1-norm of `x`.
/// * `R` is allocated:
///   * `X.has_simple_abs1` exists and is true: `fn abs(?std.mem.Allocator, X) !R`:
///     Returns the 1-norm of `x` as a newly allocated value, if the allocator
///     is provided, or a read-only view if not. If not provided, it must not
///     fail.
///   * `X.has_simple_abs1` does not exist or is false: `fn abs(std.mem.Allocator, X) !R`:
///     Returns the 1-norm of `x` as a newly allocated value.
pub inline fn abs1(x: anytype, ctx: anytype) !numeric.Abs1(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Abs1(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return int.abs(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.abs(x);
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return dyadic.abs(x);
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.abs1(x);
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
                integer.abs(ctx.allocator, x)
            else
                integer.abs(null, x) catch unreachable;
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
                rational.abs(ctx.allocator, x)
            else
                rational.abs(null, x) catch unreachable;
        },
        .real => @compileError("zml.numeric.abs1: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.abs1: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                if (comptime @hasDecl(X, "has_simple_abs1") and X.has_simple_abs1) {
                    const Impl: type = comptime types.anyHasMethod(
                        &.{ R, X },
                        "abs1",
                        fn (?std.mem.Allocator, X) anyerror!R,
                        &.{ std.mem.Allocator, X },
                    ) orelse
                        @compileError("zml.numeric.abs1: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs1(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                    return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                        Impl.abs1(ctx.allocator, x)
                    else
                        Impl.abs1(null, x) catch unreachable;
                } else {
                    const Impl: type = comptime types.anyHasMethod(
                        &.{ R, X },
                        "abs1",
                        fn (std.mem.Allocator, X) anyerror!R,
                        &.{ std.mem.Allocator, X },
                    ) orelse
                        @compileError("zml.numeric.abs1: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs1(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                    return Impl.abs1(ctx.allocator, x);
                }
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "abs1",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.abs1: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn abs1(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.abs1(x);
            }
        },
    }
}
