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

pub fn Sqrt(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.sqrt: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Sqrt(X),
        .int => return float.Sqrt(X),
        .float => return float.Sqrt(X),
        .dyadic => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return cfloat.Sqrt(X),
        .integer => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Sqrt", fn (type) type, &.{X}))
                @compileError("zml.numeric.sqrt: " ++ @typeName(X) ++ " must implement `fn Sqrt(type) type`");

            return X.Sqrt(X);
        },
    }
}

/// Returns the square root `√x` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.sqrt(x: X, ctx: anytype) !numeric.Sqrt(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the square root of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Sqrt(X)`.
///
/// #### `numeric.Sqrt(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Sqrt(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Sqrt(@TypeOf(x))`: The square root of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sqrt` method. The expected signature and
/// behavior of `Sqrt` are as follows:
/// * `fn Sqrt(type) type`: Returns the return type of `sqrt` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Sqrt(X)` as `R`. Then, `R` or `X`
/// must implement the required `sqrt` method. The expected signatures and
/// behavior of `sqrt` are as follows:
/// * `R` is not allocated: `fn sqrt(X) R`: Returns the square root of `x`.
/// * `R` is allocated: `fn sqrt(std.mem.Allocator, X) !R`: Returns the square
///   root of `x` as a newly allocated value.
pub inline fn sqrt(x: anytype, ctx: anytype) !numeric.Sqrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Sqrt(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sqrt(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sqrt(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sqrt(x);
        },
        .dyadic => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.sqrt(x);
        },
        .integer => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sqrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "sqrt",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.sqrt: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sqrt(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.sqrt(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "sqrt",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.sqrt: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sqrt(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.sqrt(x);
            }
        },
    }
}
