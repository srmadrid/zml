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

pub fn Asin(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.asin: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Asin(X),
        .int => return float.Asin(X),
        .float => return float.Asin(X),
        .dyadic => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Asin", fn (type) type, &.{X}))
                @compileError("zml.numeric.asin: " ++ @typeName(X) ++ " must implement `fn Asin(type) type`");

            return X.Asin(X);
        },
    }
}

/// Returns the arcsine `sin⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.asin(x: X, ctx: anytype) !numeric.Asin(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the arcsine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Asin(X)`.
///
/// #### `numeric.Asin(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Asin(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Asin(@TypeOf(x))`: The arcsine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Asin` method. The expected signature and
/// behavior of `Asin` are as follows:
/// * `fn Asin(type) type`: Returns the return type of `asin` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Asin(X)` as `R`. Then, `R` or `X`
/// must implement the required `asin` method. The expected signatures and
/// behavior of `asin` are as follows:
/// * `R` is not allocated: `fn asin(X) R`: Returns the arcsine of `x`.
/// * `R` is allocated: `fn asin(std.mem.Allocator, X) !R`: Returns the arcsine
///   of `x` as a newly allocated value.
pub inline fn asin(x: anytype, ctx: anytype) !numeric.Asin(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Asin(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.asin(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.asin(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.asin(x);
        },
        .dyadic => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.asin(x);
        },
        .integer => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.asin: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "asin",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.asin: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn asin(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.asin(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "asin",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.asin: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn asin(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.asin(x);
            }
        },
    }
}
