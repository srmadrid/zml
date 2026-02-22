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

pub fn Asinh(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.asinh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Asinh(X),
        .int => return float.Asinh(X),
        .float => return float.Asinh(X),
        .dyadic => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Asinh", fn (type) type, &.{X}))
                @compileError("zml.numeric.asinh: " ++ @typeName(X) ++ " must implement `fn Asinh(type) type`");

            return X.Asinh(X);
        },
    }
}

/// Returns the hyperbolic arcsine `sinh⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.asinh(x: X, ctx: anytype) !numeric.Asinh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic arcsine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Asinh(X)`.
///
/// #### `numeric.Asinh(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Asinh(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Asinh(@TypeOf(x))`: The hyperbolic arcsine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Asinh` method. The expected signature and
/// behavior of `Asinh` are as follows:
/// * `fn Asinh(type) type`: Returns the return type of `asinh` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Asinh(X)` as `R`. Then, `R` or `X`
/// must implement the required `asinh` method. The expected signatures and
/// behavior of `asinh` are as follows:
/// * `R` is not allocated: `fn asinh(X) R`: Returns the hyperbolic arcsine of
///   `x`.
/// * `R` is allocated: `fn asinh(std.mem.Allocator, X) !R`: Returns the
///   hyperbolic arcsine of `x` as a newly allocated value.
pub inline fn asinh(x: anytype, ctx: anytype) !numeric.Asinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Asinh(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.asinh(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.asinh(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.asinh(x);
        },
        .dyadic => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.asinh(x);
        },
        .integer => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.asinh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "asinh",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.asinh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn asinh(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.asinh(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "asinh",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.asinh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn asinh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.asinh(x);
            }
        },
    }
}
