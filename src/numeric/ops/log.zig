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
            if (comptime !types.hasMethod(X, "Log", fn (type) type, &.{X}))
                @compileError("zml.numeric.log: " ++ @typeName(X) ++ " must implement `fn Log(type) type`");

            return X.Log(X);
        },
    }
}

/// Returns the natural logarithm `log(x)` of a numeric `x`.
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
/// The fields of `ctx` depend on `numeric.Log(X)`.
///
/// #### `numeric.Log(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Log(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Log(@TypeOf(x))`: The natural logarithm of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Log(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Log` method. The expected signature and
/// behavior of `Log` are as follows:
/// * `fn Log(type) type`: Returns the return type of `log` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Log(X)` as `R`. Then, `R` or `X` must
/// implement the required `log` method. The expected signatures and behavior of
/// `log` are as follows:
/// * `R` is not allocated: `fn log(X) R`: Returns the logarithm of `x`.
/// * `R` is allocated: `fn log(std.mem.Allocator, X) !R`: Returns the logarithm
///   of `x` as a newly allocated value.
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
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "log",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.log: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn log(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.log(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "log",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.log: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn log(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.log(x);
            }
        },
    }
}
