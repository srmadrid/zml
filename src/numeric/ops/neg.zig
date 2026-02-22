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

pub fn Neg(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.neg: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

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
            if (comptime !types.hasMethod(X, "Neg", fn (type) type, &.{X}))
                @compileError("zml.numeric.neg: " ++ @typeName(X) ++ " must implement `fn Neg(type) type`");

            return X.Neg(X);
        },
    }
}

/// Returns the negation of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.neg(x: X, ctx: anytype) !numeric.Neg(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the negation of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Neg(X)`.
///
/// #### `numeric.Neg(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Neg(X)` is allocated
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// ## Returns
/// `numeric.Neg(@TypeOf(x))`: The negation of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Neg` method. The expected signature and
/// behavior of `Neg` are as follows:
/// * `fn Neg(type) type`: Returns the return type of `neg` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Neg(X)` as `R`. Then, `R` or `X` must
/// implement the required `neg` method. The expected signatures and behavior of
/// `neg` are as follows:
/// * `R` is not allocated: `fn neg(X) R`: Returns the negation of `x`.
/// * `R` is allocated: `fn neg(?std.mem.Allocator, X) !R`: Returns the negation
///   of `x` as a newly allocated value, if the allocator is provided, or a
///   read-only view if not. If not provided, it must not fail.
pub inline fn neg(x: anytype, ctx: anytype) !numeric.Neg(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Neg(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return !x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return -x;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return -x;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return dyadic.neg(x);
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.neg(x);
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
                integer.neg(ctx.allocator, x)
            else
                integer.neg(null, x) catch unreachable;
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
                rational.neg(ctx.allocator, x)
            else
                rational.neg(null, x) catch unreachable;
        },
        .real => @compileError("zml.numeric.neg: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.neg: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "neg",
                    fn (?std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.neg: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn neg(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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
                    Impl.neg(ctx.allocator, x)
                else
                    Impl.neg(null, x) catch unreachable;
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "neg",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.neg: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn neg(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.neg(x);
            }
        },
    }
}
