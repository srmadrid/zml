const meta = @import("../../meta.zig");

const int = @import("../../int.zig");
const float = @import("../../float.zig");
const dyadic = @import("../../dyadic.zig");
const complex = @import("../../complex.zig");

const numeric = @import("../../numeric.zig");

pub fn Re(X: type) type {
    comptime if (!meta.isNumeric(X))
        @compileError("zsl.numeric.Re: X must be a numeric type, got \n\tX = " ++ @typeName(X) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool => return X,
        .int => return X,
        .float => return X,
        .dyadic => return X,
        .complex => return meta.Scalar(X),
        .custom => {
            if (comptime !meta.hasMethod(X, "Re", fn (type) type, &.{X}))
                @compileError("zsl.numeric.Re: " ++ @typeName(X) ++ " must implement `fn Re(type) type`");

            return X.Re(X);
        },
    }
}

/// Returns the real part of a numeric `x`.
///
/// ## Signature
/// ```zig
/// numeric.re(x: X) numeric.Re(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The numeric value to get the real part of.
///
/// ## Returns
/// `numeric.Re(@TypeOf(x))`: The real part of `x`.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `X` must implement the required `Re` method. The expected signature and
/// behavior of `Re` are as follows:
/// * `fn Re(type) type`: Returns the type of the real part of `x`.
///
/// `numeric.Re(X)` or `X` must implement the required `re` method. The
/// expected signature and behavior of `re` are as follows:
/// * `fn re(X) numeric.Re(X)`: Returns the real part of `x`.
pub fn re(x: anytype) numeric.Re(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = numeric.Re(X);

    switch (comptime meta.numericType(X)) {
        .bool => return x,
        .int => return x,
        .float => return x,
        .dyadic => return x,
        .complex => return x.re,
        .custom => {
            const Impl: type = comptime meta.anyHasMethod(
                &.{ R, X },
                "re",
                fn (X) numeric.Re(X),
                &.{X},
            ) orelse
                @compileError("zsl.numeric.re: " ++ @typeName(R) ++ " or " ++ @typeName(X) ++ " must implement `fn re(" ++ @typeName(X) ++ ") " ++ @typeName(R) ++ "`");

            return Impl.re(x);
        },
    }
}
