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

pub fn Re(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.re: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return types.Scalar(X),
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return types.Scalar(X),
        .custom => {
            if (comptime !types.hasMethod(X, "Re", fn (type) type, &.{X}))
                @compileError("zml.numeric.re: " ++ @typeName(X) ++ " must implement `fn Re(type) type`");

            return X.Re(X);
        },
    }
}

/// Returns the real part of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.re(x: X, ctx: anytype) !numeric.Re(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the real part of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Re(X)`.
///
/// #### `numeric.Re(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Re(X)` is allocated
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// ## Returns
/// `numeric.Re(@TypeOf(x))`: The real part of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Re` method. The expected signature and
/// behavior of `Re` are as follows:
/// * `fn Re(type) type`: Returns the return type of `re` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Re(X)` as `R`. Then, `R` or `X` must
/// implement the required `re` method. The expected signatures and behavior of
/// `re` are as follows:
/// * `R` is not allocated: `fn re(X) R`: Returns the real part of `x`.
/// * `R` is allocated: `fn re(?std.mem.Allocator, X) !R`: Returns the real part
///   of `x` as a newly allocated value, if the allocator is provided, or a
///   read-only view if not. If not provided, it must not fail.
pub inline fn re(x: anytype, ctx: anytype) !numeric.Re(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Re(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x.re;
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
                x.copy(ctx.allocator)
            else
                x.view();
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
                x.copy(ctx.allocator)
            else
                x.view();
        },
        .real => @compileError("zml.numeric.re: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.re: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "re",
                    fn (?std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.re: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn re(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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
                    Impl.re(ctx.allocator, x)
                else
                    Impl.re(null, x) catch unreachable;
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "re",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.re: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn re(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.re(x);
            }
        },
    }
}
