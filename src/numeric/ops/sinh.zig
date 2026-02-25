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

pub fn Sinh(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.sinh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Sinh(X),
        .int => return float.Sinh(X),
        .float => return float.Sinh(X),
        .dyadic => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Sinh", fn (type) type, &.{X}))
                @compileError("zml.numeric.sinh: " ++ @typeName(X) ++ " must implement `fn Sinh(type) type`");

            return X.Sinh(X);
        },
    }
}

/// Returns the hyperbolic sine `sinh(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.sinh(x: X, ctx: anytype) !numeric.Sinh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic sine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Sinh(X)`.
///
/// #### `numeric.Sinh(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Sinh(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Sinh(@TypeOf(x))`: The hyperbolic sine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Sinh(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sinh` method. The expected signature and
/// behavior of `Sinh` are as follows:
/// * `fn Sinh(type) type`: Returns the return type of `sinh` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Sinh(X)` as `R`. Then, `R` or `X`
/// must implement the required `sinh` method. The expected signatures and
/// behavior of `sinh` are as follows:
/// * `R` is not allocated: `fn sinh(X) R`: Returns the hyperbolic sine of `x`.
/// * `R` is allocated: `fn sinh(std.mem.Allocator, X) !R`: Returns the
///   hyperbolic sine of `x` as a newly allocated value.
pub inline fn sinh(x: anytype, ctx: anytype) !numeric.Sinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Sinh(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sinh(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sinh(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.sinh(x);
        },
        .dyadic => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.sinh(x);
        },
        .integer => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.sinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "sinh",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.sinh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sinh(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.sinh(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "sinh",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.sinh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sinh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.sinh(x);
            }
        },
    }
}
