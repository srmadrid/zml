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

pub fn Exp(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.exp: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Exp(X),
        .int => return float.Exp(X),
        .float => return float.Exp(X),
        .dyadic => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Exp", fn (type) type, &.{X}))
                @compileError("zml.numeric.exp: " ++ @typeName(X) ++ " must implement `fn Exp(type) type`");

            return X.Exp(X);
        },
    }
}

/// Returns the exponential `eˣ` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.exp(x: X, ctx: anytype) !numeric.Exp(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the exponential of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Exp(X)`.
///
/// #### `numeric.Exp(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Exp(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Exp(@TypeOf(x))`: The exponential of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Exp` method. The expected signature and
/// behavior of `Exp` are as follows:
/// * `fn Exp(type) type`: Returns the return type of `exp` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Exp(X)` as `R`. Then, `R` or `X` must
/// implement the required `exp` method. The expected signatures and behavior of
/// `exp` are as follows:
/// * `R` is not allocated: `fn exp(X) R`: Returns the exponential of `x`.
/// * `R` is allocated: `fn exp(std.mem.Allocator, X) !R`: Returns the
///   exponential of `x` as a newly allocated value.
pub inline fn exp(x: anytype, ctx: anytype) !numeric.Exp(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Exp(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.exp(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.exp(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.exp(x);
        },
        .dyadic => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.exp(x);
        },
        .integer => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.exp: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "exp",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.exp: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn exp(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.exp(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "exp",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.exp: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn exp(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.exp(x);
            }
        },
    }
}
