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

pub fn Log10(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.log10: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Log10(X),
        .int => return float.Log10(X),
        .float => return float.Log10(X),
        .dyadic => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Log10", fn (type) type, &.{X}))
                @compileError("zml.numeric.log10: " ++ @typeName(X) ++ " must implement `fn Log10(type) type`");

            return X.Log10(X);
        },
    }
}

/// Returns the base-10 logarithm `log₁₀(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.log10(x: X, ctx: anytype) !numeric.Log10(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the base-10 logarithm of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `X`.
///
/// #### `X` is not allocated
/// The context must be empty.
///
/// #### `X` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Log10(@TypeOf(x))`: The base-10 logarithm of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Log10` and `log10` methods. The expected
/// signature and behavior of `Log10` and `log10` are as follows:
/// * `fn Log10(type) type`: Returns the return type of `log10` for the custom
///   numeric type.
/// * Non-allocated: `fn log10(X) X.Log10(X)`: Returns the base-10 logarithm of
///   `x`.
/// * Allocated: `fn log10(std.mem.Allocator, X) !X.Log10(X)`: Returns the
///   base-10 logarithm of `x` as a newly allocated value.
pub inline fn log10(x: anytype, ctx: anytype) !numeric.Log10(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Log10(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log10(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log10(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log10(x);
        },
        .dyadic => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.log10: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "log10", fn (std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.log10: " ++ @typeName(X) ++ " must implement `fn log10(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.log10(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "log10", fn (X) R, &.{X}))
                    @compileError("zml.numeric.log10: " ++ @typeName(X) ++ " must implement `fn log10(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.log10(x);
            }
        },
    }
}
