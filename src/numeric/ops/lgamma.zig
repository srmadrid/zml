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

pub fn Lgamma(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.lgamma: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Lgamma(X),
        .int => return float.Lgamma(X),
        .float => return float.Lgamma(X),
        .dyadic => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Lgamma", fn (type) type, &.{X}))
                @compileError("zml.numeric.lgamma: " ++ @typeName(X) ++ " must implement `fn Lgamma(type) type`");

            return X.Lgamma(X);
        },
    }
}

/// Returns the log-gamma function of a numeric `x`.
///
/// The log-gamma function is defined as:
/// $$
/// \log(\Gamma(x)) = \log\left(\int_0^\infty t^{x - 1} e^{-t} \mathrm{d}t\right).
/// $$
///
/// ## Signature
/// ```zig
/// numeric.lgamma(x: X, ctx: anytype) !numeric.Lgamma(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the log-gamma function of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Lgamma(X)`.
///
/// #### `numeric.Lgamma(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Lgamma(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Lgamma(@TypeOf(x))`: The log-gamma of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Lgamma(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Lgamma` method. The expected signature and
/// behavior of `Lgamma` are as follows:
/// * `fn Lgamma(type) type`: Returns the return type of `lgamma` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Lgamma(X)` as `R`. Then, `R` or `X`
/// must implement the required `lgamma` method. The expected signatures and
/// behavior of `lgamma` are as follows:
/// * `R` is not allocated: `fn lgamma(X) R`: Returns the log-gamma function of
///   `x`.
/// * `R` is allocated: `fn lgamma(std.mem.Allocator, X) !R`: Returns the
///   log-gamma function of `x` as a newly allocated value.
pub inline fn lgamma(x: anytype, ctx: anytype) !numeric.Lgamma(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Lgamma(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.lgamma(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.lgamma(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.lgamma(x);
        },
        .dyadic => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.lgamma: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "lgamma",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.lgamma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn lgamma(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.lgamma(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "lgamma",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.lgamma: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn lgamma(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.lgamma(x);
            }
        },
    }
}
