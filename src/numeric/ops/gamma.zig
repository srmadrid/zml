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

pub fn Gamma(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.gamma: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Gamma(X),
        .int => return float.Gamma(X),
        .float => return float.Gamma(X),
        .dyadic => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Gamma", fn (type) type, &.{X}))
                @compileError("zml.numeric.gamma: " ++ @typeName(X) ++ " must implement `fn Gamma(type) type`");

            return X.Gamma(X);
        },
    }
}

/// Returns the gamma function of a numeric `x`.
///
/// The gamma function is defined as:
/// $$
/// \Gamma(x) = \int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.gamma(x: X, ctx: anytype) !numeric.Gamma(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the gamma function of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Gamma(X)`.
///
/// #### `numeric.Gamma(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Gamma(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Gamma(@TypeOf(x))`: The gamma of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Gamma(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Gamma` method. The expected signature and
/// behavior of `Gamma` are as follows:
/// * `fn Gamma(type) type`: Returns the return type of `gamma` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Gamma(X)` as `R`. Then, `R` or `X`
/// must implement the required `gamma` method. The expected signatures and
/// behavior of `gamma` are as follows:
/// * `R` is not allocated: `fn gamma(X) R`: Returns the gamma function of `x`.
/// * `R` is allocated: `fn gamma(std.mem.Allocator, X) !R`: Returns the
///   gamma function of `x` as a newly allocated value.
pub inline fn gamma(x: anytype, ctx: anytype) !numeric.Gamma(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Gamma(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.gamma(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.gamma(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.gamma(x);
        },
        .dyadic => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.gamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "gamma",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.gamma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn gamma(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.gamma(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "gamma",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.gamma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn gamma(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.gamma(x);
            }
        },
    }
}
