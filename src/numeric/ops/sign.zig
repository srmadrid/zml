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

pub fn Sign(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.sign: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return X,
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return X,
        .custom => {
            if (comptime !types.hasMethod(X, "Sign", fn (type) type, &.{X}))
                @compileError("zml.numeric.sign: " ++ @typeName(X) ++ " must implement `fn Sign(type) type`");

            return X.Sign(X);
        },
    }
}

/// Returns the sign of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.sign(x: X, ctx: anytype) !numeric.Sign(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the sign of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Sign(X)` and `X`.
///
/// #### `numeric.Sign(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Sign(X)` is allocated and `X` is not complex
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// #### `numeric.Sign(X)` is allocated and `X` is complex
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Sign(@TypeOf(x))`: The sign of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sign` method. The expected signature and
/// behavior of `Sign` are as follows:
/// * `fn Sign(type) type`: Returns the return type of `sign` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Sign(X)` as `R`. Then, `R` or `X`
/// must implement the required `sign` method. The expected signatures and
/// behavior of `sign` are as follows:
/// * `R` is not allocated: `fn sign(X) R`: Returns the sign of `x`.
/// * `R` is allocated: `fn sign(?std.mem.Allocator, X) !R`: Returns the sign of
///   `x` as a newly allocated value, if the allocator is provided, or a
///   read-only view if not. If not provided, it must not fail.
pub inline fn sign(x: anytype, ctx: anytype) !numeric.Sign(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Sign(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return int.sign(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sign(x);
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return dyadic.sign(x);
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x.sign();
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
                integer.sign(ctx.allocator, x)
            else
                integer.sign(null, x) catch unreachable;
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
                rational.sign(ctx.allocator, x)
            else
                rational.sign(null, x) catch unreachable;
        },
        .real => @compileError("zml.numeric.sign: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sign: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                if (comptime types.isComplexType(X)) {
                    const Impl: type = comptime types.haveMethod(
                        &.{ R, X },
                        "sign",
                        fn (std.mem.Allocator, X) anyerror!R,
                        &.{ std.mem.Allocator, X },
                    ) orelse
                        @compileError("zml.numeric.sign: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sign(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                    return Impl.sign(ctx.allocator, x);
                } else {
                    const Impl: type = comptime types.haveMethod(
                        &.{ R, X },
                        "sign",
                        fn (?std.mem.Allocator, X) anyerror!R,
                        &.{ std.mem.Allocator, X },
                    ) orelse
                        @compileError("zml.numeric.sign: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sign(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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
                        Impl.sign(ctx.allocator, x)
                    else
                        Impl.sign(null, x) catch unreachable;
                }
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "sign",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.sign: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sign(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.sign(x);
            }
        },
    }
}
