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

pub fn Erfc(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.erfc: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Erfc(X),
        .int => return float.Erfc(X),
        .float => return float.Erfc(X),
        .dyadic => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Erfc", fn (type) type, &.{X}))
                @compileError("zml.numeric.erfc: " ++ @typeName(X) ++ " must implement `fn Erfc(type) type`");

            return X.Erfc(X);
        },
    }
}

/// Returns the complementary error function of a numeric `x`.
///
/// The complementary error function is defined as:
/// $$
/// \mathrm{erfc}(x) = 1 - \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.
/// $$
///
/// ## Signature
/// ```zig
/// numeric.erfc(x: X, ctx: anytype) !numeric.Erfc(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the complementary error function
///   of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Erfc(X)`.
///
/// #### `numeric.Erfc(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Erfc(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Erfc(@TypeOf(x))`: The complementary error fuction of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Erfc(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Erfc` method. The expected signature and
/// behavior of `Erfc` are as follows:
/// * `fn Erfc(type) type`: Returns the return type of `erfc` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Erfc(X)` as `R`. Then, `R` or `X`
/// must implement the required `erfc` method. The expected signatures and
/// behavior of `erfc` are as follows:
/// * `R` is not allocated: `fn erfc(X) R`: Returns the complementary error
///   function of `x`.
/// * `R` is allocated: `fn erfc(std.mem.Allocator, X) !R`: Returns the
///   complementary error function of `x` as a newly allocated value.
pub inline fn erfc(x: anytype, ctx: anytype) !numeric.Erfc(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Erfc(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.erfc(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.erfc(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.erfc(x);
        },
        .dyadic => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.erfc: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "erfc",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.erfc: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn erfc(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.erfc(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "erfc",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.erfc: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn erfc(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.erfc(x);
            }
        },
    }
}
