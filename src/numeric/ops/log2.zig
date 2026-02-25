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

pub fn Log2(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.log2: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Log2(X),
        .int => return float.Log2(X),
        .float => return float.Log2(X),
        .dyadic => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Log2", fn (type) type, &.{X}))
                @compileError("zml.numeric.log2: " ++ @typeName(X) ++ " must implement `fn Log2(type) type`");

            return X.Log2(X);
        },
    }
}

/// Returns the base-2 logarithm `log₂(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.log2(x: X, ctx: anytype) !numeric.Log2(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the base-2 logarithm of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Log2(X)`.
///
/// #### `numeric.Log2(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Log2(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Log2(@TypeOf(x))`: The base-2 logarithm of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Log2(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Log2` method. The expected signature and
/// behavior of `Log2` are as follows:
/// * `fn Log2(type) type`: Returns the return type of `log2` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Log2(X)` as `R`. Then, `R` or `X`
/// must implement the required `log2` method. The expected signatures and
/// behavior of `log2` are as follows:
/// * `R` is not allocated: `fn log2(X) R`: Returns the base-2 logarithm of `x`.
/// * `R` is allocated: `fn log2(std.mem.Allocator, X) !R`: Returns the base-2
///   logarithm of `x` as a newly allocated value.
pub inline fn log2(x: anytype, ctx: anytype) !numeric.Log2(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Log2(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log2(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log2(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log2(x);
        },
        .dyadic => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.log2: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "log2",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.log2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn log2(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.log2(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "log2",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.log2: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn log2(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.log2(x);
            }
        },
    }
}
