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

pub fn Erf(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.erf: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Erf(X),
        .int => return float.Erf(X),
        .float => return float.Erf(X),
        .dyadic => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Erf", fn (type) type, &.{X}))
                @compileError("zml.numeric.erf: " ++ @typeName(X) ++ " must implement `fn Erf(type) type`");

            return X.Erf(X);
        },
    }
}

/// Returns the error function of a numeric `x`.
///
/// The error function is defined as:
/// $$
/// \mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.erf(x: X, ctx: anytype) !numeric.Erf(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the error function of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Erf(X)`.
///
/// #### `numeric.Erf(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Erf(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Erf(@TypeOf(x))`: The error fuction of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Erf(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Erf` method. The expected signature and
/// behavior of `Erf` are as follows:
/// * `fn Erf(type) type`: Returns the return type of `erf` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Erf(X)` as `R`. Then, `R` or `X`
/// must implement the required `erf` method. The expected signatures and
/// behavior of `erf` are as follows:
/// * `R` is not allocated: `fn erf(X) R`: Returns the error function of `x`.
/// * `R` is allocated: `fn erf(std.mem.Allocator, X) !R`: Returns the error
///   function of `x` as a newly allocated value.
pub inline fn erf(x: anytype, ctx: anytype) !numeric.Erf(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Erf(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.erf(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.erf(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.erf(x);
        },
        .dyadic => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.erf: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "erf",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.erf: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn erf(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.erf(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "erf",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.erf: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn erf(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.erf(x);
            }
        },
    }
}
