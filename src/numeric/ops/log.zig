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

pub fn Log(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.log: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Log(X),
        .int => return float.Log(X),
        .float => return float.Log(X),
        .dyadic => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return cfloat.Log(X),
        .integer => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Log", fn (type) type, &.{}))
                @compileError("zml.numeric.log: " ++ @typeName(X) ++ " must implement `fn Log(type) type`");

            return X.Log(X);
        },
    }
}

/// Returns the the natural logarithm `log(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.log(x: X, ctx: anytype) !numeric.Log(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the natural logarithm of.
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
/// `numeric.Log(@TypeOf(x))`: The natural logarithm of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Log` and `log` methods. The expected
/// signature and behavior of `Log` and `log` are as follows:
/// * `fn Log(type) type`: Returns the return type of `log` for the custom
///   numeric type.
/// * Non-allocated: `fn log(X) X.Log(X)`: Returns the natural logarithm of `x`.
/// * Allocated: `fn log(std.mem.Allocator, X) !X.Log(X)`: Returns the natural
///   logarithm of `x` as a newly allocated value.
pub inline fn log(x: anytype, ctx: anytype) !numeric.Log(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Log(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.log(x);
        },
        .dyadic => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.log(x);
        },
        .integer => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.log: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "log", fn (std.mem.Allocator, X) anyerror!R, &.{}))
                    @compileError("zml.numeric.log: " ++ @typeName(X) ++ " must implement `fn log(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.log(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "log", fn (X) R, &.{}))
                    @compileError("zml.numeric.log: " ++ @typeName(X) ++ " must implement `fn log(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.log(x);
            }
        },
    }
}
