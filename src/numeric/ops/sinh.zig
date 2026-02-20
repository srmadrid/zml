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
/// The fields of `ctx` depend on `X`.
///
/// #### `X` is not allocated
/// The context must be empty.
///
/// #### `X` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Sinh(@TypeOf(x))`: The hyperbolic sine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sinh` and `sinh` methods. The expected
/// signature and behavior of `Sinh` and `sinh` are as follows:
/// * `fn Sinh(type) type`: Returns the return type of `sinh` for the custom
///   numeric type.
/// * Non-allocated: `fn sinh(X) X.Sinh(X)`: Returns the hyperbolic sine of `x`.
/// * Allocated: `fn sinh(std.mem.Allocator, X) !X.Sinh(X)`: Returns the
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
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "sinh", fn (std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.sinh: " ++ @typeName(X) ++ " must implement `fn sinh(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.sinh(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "sinh", fn (X) R, &.{X}))
                    @compileError("zml.numeric.sinh: " ++ @typeName(X) ++ " must implement `fn sinh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.sinh(x);
            }
        },
    }
}
