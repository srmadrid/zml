const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Conj(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Conj: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Conj", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Conj: " ++ @typeName(X) ++ " must implement `fn Conj(type) type`");

            return X.Conj(X);
        },
    }
}

/// Returns the complex conjugate of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.conj(x: X) numeric.Conj(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the complex conjugate of.
///
/// ## Returns
/// `numeric.Conj(@TypeOf(x))`: The complex conjugate of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Conj` method. The expected signature and
/// behavior of `Conj` are as follows:
/// * `fn Conj(type) type`: Returns the type of the complex conjugate of `x`.
///
/// `numeric.Conj(X)` or `X` must implement the required `conj` method. The
/// expected signature and behavior of `conj` are as follows:
/// * `fn conj(X) numeric.Conj(X)`: Returns the complex conjugate of `x`.
pub fn conj(x: anytype) numeric.Conj(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Conj(X);

    switch (comptime meta.numericType(X)) {
        .bool => return x,
        .int => return x,
        .float => return x,
        .dyadic => return x,
        .complex => return x.conj(),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "conj",
                fn (X) numeric.Conj(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.conj: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn conj(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.conj(x);
        },
    }
}

/// Performs computation of the complex conjugate of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.conjInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the complex conjugate of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `conjInto` method. The expected
/// signature and behavior of `conjInto` are as follows:
/// * `fn conjInto(*O, X) void`: Computes the complex conjugate of `x` and
///   stores it in `o`.
///
/// If neither `O` nor `X` implement the required `conjInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.conj`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn conjInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.conjInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "conjInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.conjInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "conjInto", fn (*O, X) void, &.{ *O, X }))
                return O.conjInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "conjInto", fn (*O, X) void, &.{ *O, X }))
            return X.conjInto(o, x);
    }

    return numeric.set(o, numeric.conj(x));
}
