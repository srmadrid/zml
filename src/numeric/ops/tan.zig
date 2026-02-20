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

pub fn Tan(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.tan: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Tan(X),
        .int => return float.Tan(X),
        .float => return float.Tan(X),
        .dyadic => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Tan", fn (type) type, &.{X}))
                @compileError("zml.numeric.tan: " ++ @typeName(X) ++ " must implement `fn Tan(type) type`");

            return X.Tan(X);
        },
    }
}

/// Returns the tangent `tan(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.tan(x: X, ctx: anytype) !numeric.Tan(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the tangent of.
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
/// `numeric.Tan(@TypeOf(x))`: The tangent of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Tan` and `tan` methods. The expected
/// signature and behavior of `Tan` and `tan` are as follows:
/// * `fn Tan(type) type`: Returns the return type of `tan` for the custom
///   numeric type.
/// * Non-allocated: `fn tan(X) X.Tan(X)`: Returns the tangent of `x`.
/// * Allocated: `fn tan(std.mem.Allocator, X) !X.Tan(X)`: Returns the tangent
///   of `x` as a newly allocated value.
pub inline fn tan(x: anytype, ctx: anytype) !numeric.Tan(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Tan(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.tan(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.tan(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.tan(x);
        },
        .dyadic => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.tan(x);
        },
        .integer => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.tan: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "tan", fn (std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.tan: " ++ @typeName(X) ++ " must implement `fn tan(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.tan(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "tan", fn (X) R, &.{X}))
                    @compileError("zml.numeric.tan: " ++ @typeName(X) ++ " must implement `fn tan(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.tan(x);
            }
        },
    }
}
