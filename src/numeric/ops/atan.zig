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

pub fn Atan(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.atan: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Atan(X),
        .int => return float.Atan(X),
        .float => return float.Atan(X),
        .dyadic => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Atan", fn (type) type, &.{X}))
                @compileError("zml.numeric.atan: " ++ @typeName(X) ++ " must implement `fn Atan(type) type`");

            return X.Atan(X);
        },
    }
}

/// Returns the arctangent `tan⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.atan(x: X, ctx: anytype) !numeric.Atan(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the arctangent of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Atan(X)`.
///
/// #### `numeric.Atan(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Atan(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Atan(@TypeOf(x))`: The arctangent of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Atan(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Atan` method. The expected signature and
/// behavior of `Atan` are as follows:
/// * `fn Atan(type) type`: Returns the return type of `atan` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Atan(X)` as `R`. Then, `R` or `X`
/// must implement the required `atan` method. The expected signatures and
/// behavior of `atan` are as follows:
/// * `R` is not allocated: `fn atan(X) R`: Returns the arctangent of `x`.
/// * `R` is allocated: `fn atan(std.mem.Allocator, X) !R`: Returns the
///   arctangent of `x` as a newly allocated value.
pub inline fn atan(x: anytype, ctx: anytype) !numeric.Atan(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Atan(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.atan(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.atan(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.atan(x);
        },
        .dyadic => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.atan(x);
        },
        .integer => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.atan: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "atan",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.atan: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atan(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.atan(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "atan",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.atan: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atan(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.atan(x);
            }
        },
    }
}
