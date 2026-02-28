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

pub fn Conj(X: type) type {
    comptime if (!types.isNumeric(X))
        @compileError("zml.numeric.conj: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime types.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .cfloat => return X,
        .integer => return X,
        .rational => return X,
        .real => return X,
        .complex => return X,
        .custom => {
            if (comptime !types.hasMethod(X, "Conj", fn (type) type, &.{X}))
                @compileError("zml.numeric.conj: " ++ @typeName(X) ++ " must implement `fn Conj(type) type`");

            return X.Conj(X);
        },
    }
}

/// Returns the complex conjugate of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.conj(x: X, ctx: anytype) !numeric.Conj(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the complex conjugate of.
/// * `ctx` (`anytype`): A context struct providing necessary resources and
///   configuration for the operation. The required fields depend on `X`. If the
///   context is missing required fields or contains unnecessary or wrongly
///   typed fields, the compiler will emit a detailed error message describing
///   the expected structure.
///
/// ### Context structure
/// The fields of `ctx` depend on `numeric.Conj(X)`.
///
/// #### `numeric.Conj(X)` is not allocated
/// The context must be empty.
///
/// #### `numeric.Conj(X)` is allocated and `X.has_simple_conj` exists and is true
/// * `allocator: std.mem.Allocator` (optional): The allocator to use for the
///   output value. If not provided, a read-only view will be returned.
///
/// #### `numeric.Conj(X)` is allocated and `X.has_simple_conj` does not exist or is false
/// * `allocator: std.mem.Allocator`: The allocator to use for the output value.
///
/// ## Returns
/// `numeric.Conj(@TypeOf(x))`: The complex conjugate of `x`.
///
/// ## Errors
/// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails. Can
///   only happen if `numeric.Conj(X)` is allocated and an allocator is
///   provided.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Conj` method. The expected signature and
/// behavior of `Conj` are as follows:
/// * `fn Conj(type) type`: Returns the return type of `conj` for the custom
///   numeric type.
///
/// Let us denote the return type `numeric.Conj(X)` as `R`. Then, `R` or `X`
/// must implement the required `conj` method. The expected signatures and
/// behavior of `conj` are as follows:
/// * `R` is not allocated: `fn conj(X) R`: Returns the conjugate of `x`.
/// * `R` is allocated: `fn conj(?std.mem.Allocator, X) !R`: Returns the
///   conjugate of `x` as a newly allocated value, if the allocator is provided,
///   or a read-only view if not. If not provided, it must not fail.
/// * `R` is allocated:
///   * `X.has_simple_conj` exists and is true: `fn conj(?std.mem.Allocator, X) !R`:
///     Returns the complex conjugate of `x` as a newly allocated value, if the
///     allocator is provided, or a read-only view if not. If not provided, it
///     must not fail.
///   * `X.has_simple_conj` does not exist or is false: `fn conj(std.mem.Allocator, X) !R`:
///     Returns the complex conjugate of `x` as a newly allocated value.
pub inline fn conj(x: anytype, ctx: anytype) !numeric.Conj(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Conj(X);

    switch (comptime types.numericType(X)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .dyadic => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x.conj();
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = false,
                        .description = "The allocator to use for the integer's memory allocation. If not provided, a view will be returned.",
                    },
                },
            );

            return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                x.copy(ctx.allocator)
            else
                x.view();
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{
                        .type = std.mem.Allocator,
                        .required = false,
                        .description = "The allocator to use for the rational's memory allocation. If not provided, a view will be returned.",
                    },
                },
            );

            return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                x.copy(ctx.allocator)
            else
                x.view();
        },
        .real => @compileError("zml.numeric.conj: not implemented for " ++ @typeName(X) ++ " yet."),
        .complex => @compileError("zml.numeric.conj: not implemented for " ++ @typeName(X) ++ " yet."),
        .custom => {
            if (comptime types.isAllocated(R)) {
                if (comptime @hasDecl(X, "has_simple_conj") and X.has_simple_conj) {
                    const Impl: type = comptime types.anyHasMethod(
                        &.{ R, X },
                        "conj",
                        fn (?std.mem.Allocator, X) anyerror!R,
                        &.{ std.mem.Allocator, X },
                    ) orelse
                        @compileError("zml.numeric.conj: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn conj(?std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{
                                .type = std.mem.Allocator,
                                .required = false,
                                .description = "The allocator to use for the custom numeric's memory allocation. If not provided, a view will be returned.",
                            },
                        },
                    );

                    return if (comptime types.ctxHasField(@TypeOf(ctx), "allocator", std.mem.Allocator))
                        Impl.conj(ctx.allocator, x)
                    else
                        Impl.conj(null, x) catch unreachable;
                } else {
                    const Impl: type = comptime types.anyHasMethod(
                        &.{ R, X },
                        "conj",
                        fn (std.mem.Allocator, X) anyerror!R,
                        &.{ std.mem.Allocator, X },
                    ) orelse
                        @compileError("zml.numeric.conj: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn conj(std.mem.Allocator, " ++ @typeName(X) ++ ") !" ++ @typeName(R) ++ "`");

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

                    return Impl.conj(ctx.allocator, x);
                }
            } else {
                const Impl: type = comptime types.anyHasMethod(
                    &.{ R, X },
                    "conj",
                    fn (X) R,
                    &.{X},
                ) orelse
                    @compileError("zml.numeric.conj: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn conj(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

                comptime types.validateContext(@TypeOf(ctx), .{});

                return Impl.conj(x);
            }
        },
    }
}
