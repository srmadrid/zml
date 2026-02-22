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

pub fn Cosh(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.cosh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Cosh(X),
        .int => return float.Cosh(X),
        .float => return float.Cosh(X),
        .dyadic => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Cosh", fn (type) type, &.{X}))
                @compileError("zml.numeric.cosh: " ++ @typeName(X) ++ " must implement `fn Cosh(type) type`");

            return X.Cosh(X);
        },
    }
}

/// Returns the hyperbolic cosine `cosh(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.cosh(x: X, ctx: anytype) !numeric.Cosh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic cosine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Cosh(X)`.
///
/// #### `numeric.Cosh(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Cosh(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Cosh(@TypeOf(x))`: The hyperbolic cosine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Cosh` method. The expected signature and
/// behavior of `Cosh` are as follows:
/// * `fn Cosh(type) type`: Returns the return type of `cosh` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Cosh(X)` as `R`. Then, `R` or `X`
/// must implement the required `cosh` method. The expected signatures and
/// behavior of `cosh` are as follows:
/// * `R` is not allocated: `fn cosh(X) R`: Returns the hyperbolic cosine of
///   `x`.
/// * `R` is allocated: `fn cosh(std.mem.Allocator, X) !R`: Returns the
///   hyperbolic cosine of `x` as a newly allocated value.
pub inline fn cosh(x: anytype, ctx: anytype) !numeric.Cosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Cosh(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cosh(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cosh(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cosh(x);
        },
        .dyadic => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.cosh(x);
        },
        .integer => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.cosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "cosh",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.cosh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cosh(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.cosh(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "cosh",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.cosh: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cosh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.cosh(x);
            }
        },
    }
}
