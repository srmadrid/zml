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

pub fn Cbrt(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.cbrt: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Cbrt(X),
        .int => return float.Cbrt(X),
        .float => return float.Cbrt(X),
        .dyadic => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Cbrt", fn (type) type, &.{X}))
                @compileError("zml.numeric.cbrt: " ++ @typeName(X) ++ " must implement `fn Cbrt(type) type`");

            return X.Cbrt(X);
        },
    }
}

/// Returns the cube root `∛x` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.cbrt(x: X, ctx: anytype) !numeric.Cbrt(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the cube root of.
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
/// `numeric.Cbrt(@TypeOf(x))`: The cube root of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Cbrt` and `cbrt` methods. The expected
/// signature and behavior of `Cbrt` and `cbrt` are as follows:
/// * `fn Cbrt(type) type`: Returns the return type of `cbrt` for the custom
///   numeric type.
/// * Non-allocated: `fn cbrt(X) X.Cbrt(X)`: Returns the cube root of `x`.
/// * Allocated: `fn cbrt(std.mem.Allocator, X) !X.Cbrt(X)`: Returns the cube
///   root of `x` as a newly allocated value.
pub inline fn cbrt(x: anytype, ctx: anytype) !numeric.Cbrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Cbrt(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cbrt(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cbrt(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.cbrt(x);
        },
        .dyadic => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .integer => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.cbrt: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "cbrt", fn (std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.cbrt: " ++ @typeName(X) ++ " must implement `fn cbrt(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.cbrt(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "cbrt", fn (X) R, &.{X}))
                    @compileError("zml.numeric.cbrt: " ++ @typeName(X) ++ " must implement `fn cbrt(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.cbrt(x);
            }
        },
    }
}
