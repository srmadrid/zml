const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Atan(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.atan: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.atan: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.atan: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Atan", fn (type) type, &.{X}))
                @compileError("zsl.numeric.atan: " ++ @typeName(X) ++ " must implement `fn Atan(type) type`");

            return X.Atan(X);
        },
    }
}

/// Returns the arctangent `tan⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.atan(x: X) numeric.Atan(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the arctangent of.
///
/// ## Returns
/// `numeric.Atan(@TypeOf(x))`: The arctangent of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Atan` method. The expected signature and
/// behavior of `Atan` are as follows:
/// * `fn Atan(type) type`: Returns the type of the arctangent of `x`.
///
/// `numeric.Atan(X)` or `X` must implement the required `atan` method. The
/// expected signature and behavior of `atan` are as follows:
/// * `fn atan(X) numeric.Atan(X)`: Returns the arctangent of `x`.
pub fn atan(x: anytype) numeric.Atan(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Atan(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.atan(x),
        .dyadic => return dyadic.atan(x),
        .complex => return complex.atan(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "atan",
                fn (X) numeric.Atan(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.atan: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn atan(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.atan(x);
        },
    }
}

/// Performs computation of the arctangent `tan⁻¹(x)` of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.atanInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the arctangent of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `atanInto` method. The expected
/// signature and behavior of `atanInto` are as follows:
/// * `fn atanInto(*O, X) void`: Computes the arctangent of `x` and stores it
///   in `o`.
///
/// If neither `O` nor `X` implement the required `atanInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.atan`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn atanInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.atanInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "atanInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.atanInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "atanInto", fn (*O, X) void, &.{ *O, X }))
                return O.atanInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "atanInto", fn (*O, X) void, &.{ *O, X }))
            return X.atanInto(o, x);
    }

    return numeric.set(o, numeric.atan(x));
}
