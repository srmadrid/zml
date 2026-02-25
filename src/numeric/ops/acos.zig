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

pub fn Acos(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.acos: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return float.Acos(X),
        .int => return float.Acos(X),
        .float => return float.Acos(X),
        .dyadic => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => return X,
        .integer => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime !types.hasMethod(X, "Acos", fn (type) type, &.{X}))
                @compileError("zml.numeric.acos: " ++ @typeName(X) ++ " must implement `fn Acos(type) type`");

            return X.Acos(X);
        },
    }
}

/// Returns the arccosine `cos⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.acos(x: X, ctx: anytype) !numeric.Acos(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the arccosine of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Acos(X)`.
///
/// #### `numeric.Acos(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Acos(X)` is allocated
/// * `allocator: std.mem.Allocator` The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Acos(@TypeOf(x))`: The arccosine of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Acos(X)` is allocated.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Acos` method. The expected signature and
/// behavior of `Acos` are as follows:
/// * `fn Acos(type) type`: Returns the return type of `acos` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Acos(X)` as `R`. Then, `R` or `X`
/// must implement the required `acos` method. The expected signatures and
/// behavior of `acos` are as follows:
/// * `R` is not allocated: `fn acos(X) R`: Returns the arccosine of `x`.
/// * `R` is allocated: `fn acos(std.mem.Allocator, X) !R`: Returns the
///   arccosine of `x` as a newly allocated value.
pub inline fn acos(x: anytype, ctx: anytype) !numeric.Acos(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Acos(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.acos(x);
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.acos(x);
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return float.acos(x);
        },
        .dyadic => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return cfloat.acos(x);
        },
        .integer => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .rational => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .real => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.acos: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "acos",
                    fn (std.mem.Allocator, X) anyerror!R,
                    &.{ std.mem.Allocator, X },
                ) orelse
                    @compileError("zml.numeric.acos: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn acos(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                return Impl.acos(ctx.allocator, x);
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "acos",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.acos: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn acos(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.acos(x);
            }
        },
    }
}
