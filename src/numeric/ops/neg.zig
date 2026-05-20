const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Neg(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.neg: x must be a numeric, got \n\tx: " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return X,
        .custom => {
            if (comptime !meta.hasMethod(X, "Neg", fn (type) type, &.{X}))
                @compileError("zsl.numeric.neg: " ++ @typeName(X) ++ " must implement `fn Neg(type) type`");

            return X.Neg(X);
        },
    }
}

/// Returns the negation of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.neg(x: X) numeric.Neg(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the negation of.
///
/// ## Returns
/// `numeric.Neg(@TypeOf(x))`: The negation of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Neg` method. The expected signature and
/// behavior of `Neg` are as follows:
/// * `fn Neg(type) type`: Returns the type of the negation of `x`.
///
/// `numeric.Neg(X)` or `X` must implement the required `neg` method. The
/// expected signature and behavior of `neg` are as follows:
/// * `fn neg(X) numeric.Neg(X)`: Returns the negation of `x`.
pub fn neg(x: anytype) numeric.Neg(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Neg(X);

    switch (comptime meta.numericType(X)) {
        .bool => return !x,
        .int => return -x,
        .float => return -x,
        .dyadic => return x.neg(),
        .complex => return x.neg(),
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "neg",
                fn (X) numeric.Neg(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.neg: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn neg(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.neg(x);
        },
    }
}
