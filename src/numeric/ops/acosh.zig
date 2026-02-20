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

pub fn Acosh(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.acosh: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Acosh(X),
        .int => return float.Acosh(X),
        .float => return float.Acosh(X),
        .dyadic => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Acosh", fn (type) type, &.{X}))
                @compileError("zml.numeric.acosh: " ++ @typeName(X) ++ " must implement `fn Acosh(type) type`");

            return X.Acosh(X);
        },
    }
}

/// Returns the hyperbolic arccosine `cosh⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.acosh(x: X, ctx: anytype) !numeric.Acosh(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the hyperbolic arccosine of.
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
/// `numeric.Acosh(@TypeOf(x))`: The hyperbolic arccosine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `X` is allocated and an allocator is provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Acosh` and `acosh` methods. The expected
/// signature and behavior of `Acosh` and `acosh` are as follows:
/// * `fn Acosh(type) type`: Returns the return type of `acosh` for the custom
///   numeric type.
/// * Non-allocated: `fn acosh(X) X.Acosh(X)`: Returns the hyperbolic arccosine
///   of `x`.
/// * Allocated: `fn acosh(std.mem.Allocator, X) !X.Acosh(X)`: Returns the
///   hyperbolic arccosine of `x` as a newly allocated value.
pub inline fn acosh(x: anytype, ctx: anytype) !numeric.Acosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Acosh(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.acosh(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.acosh(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.acosh(x);
        },
        .dyadic => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.acosh(x);
        },
        .integer => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.acosh: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(X)) {
                comptime if (!types.hasMethod(X, "acosh", fn (std.mem.Allocator, X) anyerror!R, &.{ std.mem.Allocator, X }))
                    @compileError("zml.numeric.acosh: " ++ @typeName(X) ++ " must implement `fn acosh(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return X.acosh(ctx.allocator, x);
            } else {
                comptime if (!types.hasMethod(X, "acosh", fn (X) R, &.{X}))
                    @compileError("zml.numeric.acosh: " ++ @typeName(X) ++ " must implement `fn acosh(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return X.acosh(x);
            }
        },
    }
}
