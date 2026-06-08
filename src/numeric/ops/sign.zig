const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Sign(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Sign: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Sign", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Sign: " ++ @typeName(X) ++ " must implement `fn Sign(type) type`");

            return X.Sign(X);
        },
    }
}

/// Returns the sign of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.sign(x: X) numeric.Sign(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the sign of.
///
/// ## Returns
/// `numeric.Sign(@TypeOf(x))`: The sign of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Sign` method. The expected signature and
/// behavior of `Sign` are as follows:
/// * `fn Sign(type) type`: Returns the type of the sign of `x`.
///
/// `numeric.Sign(X)` or `X` must implement the required `sign` method. The
/// expected signature and behavior of `sign` are as follows:
/// * `fn sign(X) numeric.Sign(X)`: Returns the sign of `x`.
pub fn sign(x: anytype) numeric.Sign(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Sign(X);

    switch (comptime meta.numericType(X)) {
        .bool => return x,
        .int => return int.sign(x),
        .float => return float.sign(x),
        .dyadic => return dyadic.sign(x),
        .complex => return complex.sign(x),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "sign",
                fn (X) numeric.Sign(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.sign: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn sign(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.sign(x);
        },
    }
}
