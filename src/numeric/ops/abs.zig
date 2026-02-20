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

pub fn Abs(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.abs: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return cfloat.Abs(X),
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return complex.Abs(X),
        .custom => {
            if (comptime !types.hasMethod(X, "Abs", fn (type) type, &.{X}))
                @compileError("zml.numeric.abs: " ++ @typeName(X) ++ " must implement `fn Abs(type) type`");

            return X.Abs(X);
        },
    }
}

/// Returns the absolute value of a numeric `x`. For complex values, returns the
/// magnitude.
///
/// ## Signature
/// ```zig
/// numeric.abs(x: X, ctx: anytype) !numeric.Abs(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the absolute value of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `X`.
///
/// #### `X` is not allocated
/// The context must be empty.
///
/// #### `X` is allocated
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// ## Returns
/// `numeric.Abs(@TypeOf(x))`: The absolute value of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Abs` and `abs` methods. The expected
/// signature and behavior of `Abs` and `abs` are as follows:
/// * `fn Abs(type) type`: Returns the return type of `abs` for the custom
///   numeric type.
/// * Non-allocated: `fn abs(X) X.Abs(X)`: Returns the absolute value of `x`.
/// * Allocated: `fn abs(?std.mem.Allocator, X) !X.Abs(X)`: Returns the absolute
///   value of `x` as a newly allocated value, if the allocator is provided, or
///   a read-only view if not.
pub inline fn abs(x: anytype, ctx: anytype) !numeric.Abs(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Abs(X);

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

            return cfloat.abs(x);
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
        .real => @compileError("zml.numeric.abs: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.abs: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "abs", fn (?std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.abs: " ++ @typeName(X) ++ " must implement `fn abs(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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
                    X.abs(ctx.allocator, x)
                else
                    X.abs(null, x);
            } else {
                comptime if (!types.hasMethod(X, "abs", fn (X) R, &.{X}))
                    @compileError("zml.numeric.abs: " ++ @typeName(X) ++ " must implement `fn abs(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.abs(x);
            }
        },
    }
}
