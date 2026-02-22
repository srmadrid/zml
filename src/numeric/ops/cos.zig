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

pub fn Cos(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.cos: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Cos(X),
        .int => return float.Cos(X),
        .float => return float.Cos(X),
        .dyadic => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Cos", fn (type) type, &.{X}))
                @compileError("zml.numeric.cos: " ++ @typeName(X) ++ " must implement `fn Cos(type) type`");

            return X.Cos(X);
        },
    }
}

/// Returns the cosine `cos(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.cos(x: X, ctx: anytype) !numeric.Cos(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the cosine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Cos(X)`.
///
/// #### `numeric.Cos(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Cos(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Cos(@TypeOf(x))`: The cosine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Cos` method. The expected signature and
/// behavior of `Cos` are as follows:
/// * `fn Cos(type) type`: Returns the return type of `cos` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Cos(X)` as `R`. Then, `R` or `X`
/// must implement the required `cos` method. The expected signatures and
/// behavior of `cos` are as follows:
/// * `R` is not allocated: `fn cos(X) R`: Returns the cosine of `x`.
/// * `R` is allocated: `fn cos(std.mem.Allocator, X) !R`: Returns the cosine of
///   `x` as a newly allocated value.
pub inline fn cos(x: anytype, ctx: anytype) !numeric.Cos(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Cos(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cos(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cos(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cos(x);
        },
        .dyadic => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.cos(x);
        },
        .integer => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.cos: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "cos",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.cos: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cos(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.cos(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.haveMethod(
                    &.{ R, X },
                    "cos",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.cos: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn cos(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.cos(x);
            }
        },
    }
}
