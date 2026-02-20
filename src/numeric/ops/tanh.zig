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

pub fn Tanh(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.tanh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Tanh(X),
        .int => return float.Tanh(X),
        .float => return float.Tanh(X),
        .dyadic => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Tanh", fn (type) type, &.{X}))
                @compileError("zml.numeric.tanh: " ++ @typeName(X) ++ " must implement `fn Tanh(type) type`");

            return X.Tanh(X);
        },
    }
}

/// Returns the hyperbolic tangent `tanh(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.tanh(x: X, ctx: anytype) !numeric.Tanh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic tangent of.
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
/// `numeric.Tanh(@TypeOf(x))`: The hyperbolic tangent of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Tanh` and `tanh` methods. The expected
/// signature and behavior of `Tanh` and `tanh` are as follows:
/// * `fn Tanh(type) type`: Returns the return type of `tanh` for the custom
///   numeric type.
/// * Non-allocated: `fn tanh(X) X.Tanh(X)`: Returns the hyperbolic tangent of
///   `x`.
/// * Allocated: `fn tanh(std.mem.Allocator, X) !X.Tanh(X)`: Returns the
///   hyperbolic tangent of `x` as a newly allocated value.
pub inline fn tanh(x: anytype, ctx: anytype) !numeric.Tanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Tanh(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.tanh(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.tanh(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.tanh(x);
        },
        .dyadic => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.tanh(x);
        },
        .integer => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.tanh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "tanh", fn (std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.tanh: " ++ @typeName(X) ++ " must implement `fn tanh(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.tanh(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "tanh", fn (X) R, &.{X}))
                    @compileError("zml.numeric.tanh: " ++ @typeName(X) ++ " must implement `fn tanh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.tanh(x);
            }
        },
    }
}
