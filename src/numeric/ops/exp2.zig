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

pub fn Exp2(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.exp2: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Exp2(X),
        .int => return float.Exp2(X),
        .float => return float.Exp2(X),
        .dyadic => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Exp2", fn (type) type, &.{X}))
                @compileError("zml.numeric.exp2: " ++ @typeName(X) ++ " must implement `fn Exp2(type) type`");

            return X.Exp2(X);
        },
    }
}

/// Returns the base-2 exponential `2ˣ` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.exp2(x: X, ctx: anytype) !numeric.Exp2(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the base-2 exponential of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Exp(X)`.
///
/// #### `numeric.Exp2(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Exp2(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Exp2(@TypeOf(x))`: The base-2 exponential of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Exp2` method. The expected signature and
/// behavior of `Exp2` are as follows:
/// * `fn Exp2(type) type`: Returns the return type of `exp2` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Exp2(X)` as `R`. Then, `R` or `X` must
/// implement the required `exp2` method. The expected signatures and behavior of
/// `exp2` are as follows:
/// * `R` is not allocated: `fn exp2(X) R`: Returns the base-2 exponential of
///   `x`.
/// * `R` is allocated: `fn exp2(std.mem.Allocator, X) !R`: Returns the base-2
///   exponential of `x` as a newly allocated value.
pub inline fn exp2(x: anytype, ctx: anytype) !numeric.Exp2(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Exp2(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.exp2(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.exp2(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.exp2(x);
        },
        .dyadic => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.exp2: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "exp2",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.exp2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn exp2(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.exp2(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "exp2",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.exp2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn exp2(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.exp2(x);
            }
        },
    }
}
