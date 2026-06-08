const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Acos(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Acos: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => @compileError("zsl.numeric.Acos: not defined for " ++ @typeName(X) ++ "."),
        .int => @compileError("zsl.numeric.Acos: not defined for " ++ @typeName(X) ++ "."),
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Acos", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Acos: " ++ @typeName(X) ++ " must implement `fn Acos(type) type`");

            return X.Acos(X);
        },
    }
}

/// Returns the arccosine `cos⁻¹(x)` of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.acos(x: X) numeric.Acos(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the arccosine of.
///
/// ## Returns
/// `numeric.Acos(@TypeOf(x))`: The arccosine of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Acos` method. The expected signature and
/// behavior of `Acos` are as follows:
/// * `fn Acos(type) type`: Returns the type of the arccosine of `x`.
///
/// `numeric.Acos(X)` or `X` must implement the required `acos` method. The
/// expected signature and behavior of `acos` are as follows:
/// * `fn acos(X) numeric.Acos(X)`: Returns the arccosine of `x`.
pub fn acos(x: anytype) numeric.Acos(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Acos(X);

    switch (comptime meta.numericType(X)) {
        .bool => unreachable,
        .int => unreachable,
        .float => return float.acos(x),
        .dyadic => return dyadic.acos(x),
        .complex => return complex.acos(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "acos",
                fn (X) numeric.Acos(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.acos: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn acos(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.acos(x);
        },
    }
}

/// Performs computation of the arccosine `cos⁻¹(x)` of a numeric `x` into a
/// numeric `o`.
///
/// ## Signature
/// ```zig
/// numeric.acosInto(o: *O, x: X) void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The numeric value to get the arccosine of.
///
/// ## Returns
/// `void`
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `O` or `X` should implement the required `acosInto` method. The expected
/// signature and behavior of `acosInto` are as follows:
/// * `fn acosInto(*O, X) void`: Computes the arccosine of `x` and stores it in
///   `o`.
///
/// If neither `O` nor `X` implement the required `acosInto` method, the
/// function will fall back to using `numeric.set` with the result of
/// `numeric.acos`, potentially resulting in a less efficient implementation. In
/// this case, `O` and `X` must adhere to the requirements of these functions.
pub fn acosInto(o: anytype, x: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or
        !meta.isNumeric(meta.Child(O)) or
        !meta.isNumeric(X))
        @compileError("zsl.numeric.acosInto: o must be a mutable one-item pointer to a numeric, and x must be a numeric, got \n\to: " ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    if (comptime meta.isCustomNumeric(O)) {
        if (comptime meta.isCustomNumeric(X)) { // O and X both custom
            if (comptime meta.anyHasMethod(&.{ O, X }, "acosInto", fn (*O, X) void, &.{ *O, X })) |Impl|
                return Impl.acosInto(o, x);
        } else { // only O custom
            if (comptime meta.hasMethod(O, "acosInto", fn (*O, X) void, &.{ *O, X }))
                return O.acosInto(o, x);
        }
    } else if (comptime meta.isCustomNumeric(X)) { // only X custom
        if (comptime meta.hasMethod(X, "acosInto", fn (*O, X) void, &.{ *O, X }))
            return X.acosInto(o, x);
    }

    return numeric.set(o, numeric.acos(x));
}
