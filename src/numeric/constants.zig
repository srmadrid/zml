const meta = @import("../meta.zig");

const int = @import("../int.zig");
const float = @import("../float.zig");
const dyadic = @import("../dyadic.zig");
const complex = @import("../complex.zig");

/// Returns the maximum representable finite value (highest point on the number
/// line) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the highest value for.
///
/// ## Returns
/// `N`: The maximum representable value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `highest: N` declaration, or implement the `highest`
/// method. The expected signature and behavior of `highest` are as follows:
/// * `fn highest(anytype) N`: Returns the highest representable value.
pub fn highest(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.highest: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return int.highest(N),
        .float => return float.highest(N),
        .dyadic => return dyadic.highest(N),
        .complex => @compileError("zsl.numeric.highest: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "highest") and @TypeOf(N.highest) == N)
                return N.highest
            else if (comptime meta.hasMethod(N, "highest", fn () N, &.{}))
                return N.highest()
            else
                @compileError("zsl.numeric.highest: " ++ @typeName(N) ++ " must expose a `highest: " ++ @typeName(N) ++ "` declaration, or implement `fn highest() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the minimum representable finite value (lowest point on the number
/// line) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the lowest value for.
///
/// ## Returns
/// `N`: The minimum representable value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `lowest: N` declaration, or implement the `lowest`
/// method. The expected signature and behavior of `lowest` are as follows:
/// * `fn lowest(anytype) N`: Returns the lowest representable value.
pub fn lowest(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.lowest: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return int.lowest(N),
        .float => return float.lowest(N),
        .dyadic => return dyadic.lowest(N),
        .complex => @compileError("zsl.numeric.lowest: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "lowest") and @TypeOf(N.lowest) == N)
                return N.lowest
            else if (comptime meta.hasMethod(N, "lowest", fn () N, &.{}))
                return N.lowest()
            else
                @compileError("zsl.numeric.lowest: " ++ @typeName(N) ++ " must expose a `lowest: " ++ @typeName(N) ++ "` declaration, or implement `fn lowest() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the smallest positive magnitude strictly greater than zero (closest
/// to zero) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the smallest positive value
///   for.
///
/// ## Returns
/// `N`: The smallest positive non-zero magnitude.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `smallest: N` declaration, or implement the `smallest`
/// method. The expected signature and behavior of `smallest` are as follows:
/// * `fn smallest(anytype) N`: Returns the smallest positive non-zero magnitude.
pub fn smallest(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.smallest: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return true,
        .int => return int.smallest(N),
        .float => return float.smallest(N),
        .dyadic => return dyadic.smallest(N),
        .complex => @compileError("zsl.numeric.smallest: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "smallest") and @TypeOf(N.smallest) == N)
                return N.smallest
            else if (comptime meta.hasMethod(N, "smallest", fn () N, &.{}))
                return N.smallest()
            else
                @compileError("zsl.numeric.smallest: " ++ @typeName(N) ++ " must expose a `smallest: " ++ @typeName(N) ++ "` declaration, or implement `fn smallest() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the machine epsilon (`ε`) for the given numeric type `N`. This
/// represents the upper bound on the relative approximation error due to
/// rounding, typically defined as the difference between `1` and the next
/// representable value (i.e., the smallest positive value `ε > 0` such that
/// `1 + ε ≠ 1`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the machine epsilon for.
///
/// ## Returns
/// `N`: The machine epsilon (`ε`).
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `eps: N` declaration, or implement the `eps`
/// method. The expected signature and behavior of `eps` are as follows:
/// * `fn eps(anytype) N`: Returns the machine epsilon value.
pub fn eps(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.eps: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => return false,
        .int => return int.eps(N),
        .float => return float.eps(N),
        .dyadic => return dyadic.eps(N),
        .complex => return complex.eps(N), // ε + 0i
        .custom => {
            if (comptime @hasDecl(N, "eps") and @TypeOf(N.eps) == N)
                return N.eps
            else if (comptime meta.hasMethod(N, "eps", fn () N, &.{}))
                return N.eps()
            else
                @compileError("zsl.numeric.eps: " ++ @typeName(N) ++ " must expose a `eps: " ++ @typeName(N) ++ "` declaration, or implement `fn eps() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns positive infinity (`∞`) for the given numeric type `N`.
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the infinity value for.
///
/// ## Returns
/// `N`: The positive infinity value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `inf: N` declaration, or implement the `inf`
/// method. The expected signature and behavior of `inf` are as follows:
/// * `fn inf(anytype) N`: Returns the positive infinity value.
pub fn inf(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.inf: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.inf: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.inf: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.inf(N),
        .dyadic => return dyadic.inf(N),
        .complex => @compileError("zsl.numeric.inf: not defined for " ++ @typeName(N) ++ "."),
        .custom => {
            if (comptime @hasDecl(N, "inf") and @TypeOf(N.inf) == N)
                return N.inf
            else if (comptime meta.hasMethod(N, "inf", fn () N, &.{}))
                return N.inf()
            else
                @compileError("zsl.numeric.inf: " ++ @typeName(N) ++ " must expose a `inf: " ++ @typeName(N) ++ "` declaration, or implement `fn inf() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the Not-a-Number (NaN) value for the given numeric type `N`. NaN
/// is used to represent undefined or unrepresentable results in mathematical
/// operations (e.g., 0 ÷ 0).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the NaN value for.
///
/// ## Returns
/// `N`: The Not-a-Number value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `nan: N` declaration, or implement the `nan`
/// method. The expected signature and behavior of `nan` are as follows:
/// * `fn nan(anytype) N`: Returns the Not-a-Number value.
pub fn nan(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.nan: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.nan: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.nan: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.nan(N),
        .dyadic => return dyadic.nan(N),
        .complex => return complex.nan(N), // nan + nan i
        .custom => {
            if (comptime @hasDecl(N, "nan") and @TypeOf(N.nan) == N)
                return N.nan
            else if (comptime meta.hasMethod(N, "nan", fn () N, &.{}))
                return N.nan()
            else
                @compileError("zsl.numeric.nan: " ++ @typeName(N) ++ " must expose a `nan: " ++ @typeName(N) ++ "` declaration, or implement `fn nan() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the mathematical constant pi (`π`) for the given numeric type `N`.
/// This represents the ratio of a circle's circumference to its diameter
/// (`C/d ≈ 3.14159…`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `π` vaalue for.
///
/// ## Returns
/// `N`: The `π` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `pi: N` declaration, or implement the `pi`
/// method. The expected signature and behavior of `pi` are as follows:
/// * `fn pi(anytype) N`: Returns the `π` value.
pub fn pi(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.pi: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.pi: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.pi: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.pi(N),
        .dyadic => return dyadic.pi(N),
        .complex => return complex.pi(N), // pi + 0i
        .custom => {
            if (comptime @hasDecl(N, "pi") and @TypeOf(N.pi) == N)
                return N.pi
            else if (comptime meta.hasMethod(N, "pi", fn () N, &.{}))
                return N.pi()
            else
                @compileError("zsl.numeric.pi: " ++ @typeName(N) ++ " must expose a `pi: " ++ @typeName(N) ++ "` declaration, or implement `fn pi() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the mathematical constant tau (`τ`) for the given numeric type `N`.
/// This represents the ratio of a circle's circumference to its radius
/// (`C/r ≈ 6.28318…`, equivalent to `2π`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `τ` value for.
///
/// ## Returns
/// `N`: The `τ` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `tau: N` declaration, or implement the `tau`
/// method. The expected signature and behavior of `tau` are as follows:
/// * `fn tau(anytype) N`: Returns the `τ` value.
pub fn tau(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.tau: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.tau: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.tau: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.tau(N),
        .dyadic => return dyadic.tau(N),
        .complex => return complex.tau(N), // tau + 0i
        .custom => {
            if (comptime @hasDecl(N, "tau") and @TypeOf(N.tau) == N)
                return N.tau
            else if (comptime meta.hasMethod(N, "tau", fn () N, &.{}))
                return N.tau()
            else
                @compileError("zsl.numeric.tau: " ++ @typeName(N) ++ " must expose a `tau: " ++ @typeName(N) ++ "` declaration, or implement `fn tau() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns Euler's number (`e`) for the given numeric type `N`. This represents
/// the base of the natural logarithm (`∑ₖ₌₀^∞ 1/k! ≈ 2.71828…`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `e` value for.
///
/// ## Returns
/// `N`: The `e` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `e: N` declaration, or implement the `e`
/// method. The expected signature and behavior of `e` are as follows:
/// * `fn e(anytype) N`: Returns the `e` value.
pub fn e(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.e: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.e: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.e: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.e(N),
        .dyadic => return dyadic.e(N),
        .complex => return complex.e(N), // e + 0i
        .custom => {
            if (comptime @hasDecl(N, "e") and @TypeOf(N.e) == N)
                return N.e
            else if (comptime meta.hasMethod(N, "e", fn () N, &.{}))
                return N.e()
            else
                @compileError("zsl.numeric.e: " ++ @typeName(N) ++ " must expose a `e: " ++ @typeName(N) ++ "` declaration, or implement `fn e() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the golden ratio (`φ`) for the given numeric type `N`. This
/// represents the positive solution to the equation `x² - x - 1 = 0`
/// (`(1 + √5) / 2 ≈ 1.61803…`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `φ` value for.
///
/// ## Returns
/// `N`: The `φ` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `phi: N` declaration, or implement the `phi`
/// method. The expected signature and behavior of `phi` are as follows:
/// * `fn phi(anytype) N`: Returns the `φ` value.
pub fn phi(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.phi: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.phi: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.phi: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.phi(N),
        .dyadic => return dyadic.phi(N),
        .complex => return complex.phi(N), // phi + 0i
        .custom => {
            if (comptime @hasDecl(N, "phi") and @TypeOf(N.phi) == N)
                return N.phi
            else if (comptime meta.hasMethod(N, "phi", fn () N, &.{}))
                return N.phi()
            else
                @compileError("zsl.numeric.phi: " ++ @typeName(N) ++ " must expose a `phi: " ++ @typeName(N) ++ "` declaration, or implement `fn phi() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns the Euler-Mascheroni constant (`γ`) for the given numeric type `N`.
/// This represents the limiting difference between the harmonic series and the
/// natural logarithm (`limₙ→∞ (∑ₖ₌₁ⁿ 1/k - ln(n)) ≈ 0.57721…`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `γ` value for.
///
/// ## Returns
/// `N`: The `γ` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `egamma: N` declaration, or implement the `egamma`
/// method. The expected signature and behavior of `egamma` are as follows:
/// * `fn egamma(anytype) N`: Returns the `γ` value.
pub fn egamma(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.egamma: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.egamma: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.egamma: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.egamma(N),
        .dyadic => return dyadic.egamma(N),
        .complex => return complex.egamma(N), // egamma + 0i
        .custom => {
            if (comptime @hasDecl(N, "egamma") and @TypeOf(N.egamma) == N)
                return N.egamma
            else if (comptime meta.hasMethod(N, "egamma", fn () N, &.{}))
                return N.egamma()
            else
                @compileError("zsl.numeric.egamma: " ++ @typeName(N) ++ " must expose a `egamma: " ++ @typeName(N) ++ "` declaration, or implement `fn egamma() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns Catalan's constant (`G`) for the given numeric type `N`. This
/// represents the alternating sum of the reciprocals of the odd square numbers
/// (`∑ₙ₌₀^∞ (-1)ⁿ/(2n + 1)² ≈ 0.91596…`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `G` value for.
///
/// ## Returns
/// `N`: The `G` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `catalan: N` declaration, or implement the `catalan`
/// method. The expected signature and behavior of `catalan` are as follows:
/// * `fn catalan(anytype) N`: Returns the `G` value.
pub fn catalan(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.catalan: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.catalan: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.catalan: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.catalan(N),
        .dyadic => return dyadic.catalan(N),
        .complex => return complex.catalan(N), // catalan + 0i
        .custom => {
            if (comptime @hasDecl(N, "catalan") and @TypeOf(N.catalan) == N)
                return N.catalan
            else if (comptime meta.hasMethod(N, "catalan", fn () N, &.{}))
                return N.catalan()
            else
                @compileError("zsl.numeric.catalan: " ++ @typeName(N) ++ " must expose a `catalan: " ++ @typeName(N) ++ "` declaration, or implement `fn catalan() " ++ @typeName(N) ++ "`");
        },
    }
}

/// Returns Apéry's constant (`ζ(3)`) for the given numeric type `N`. This
/// represents the value of the Riemann zeta function evaluated at 3
/// (`∑ₖ₌₁^∞ 1/k³ ≈ 1.202056…`).
///
/// ## Arguments
/// * `N` (`comptime type`): The type to generate the `ζ(3)` value for.
///
/// ## Returns
/// `N`: The `ζ(3)` value.
///
/// ## Custom type support
/// This function supports custom numeric types via specific method
/// implementations.
///
/// `N` must expose the `apery: N` declaration, or implement the `apery`
/// method. The expected signature and behavior of `apery` are as follows:
/// * `fn apery(anytype) N`: Returns the `ζ(3)` value.
pub fn apery(comptime N: type) N {
    comptime if (!meta.isNumeric(N))
        @compileError("zsl.numeric.apery: " ++ @typeName(N) ++ " is not a numeric type");

    switch (comptime meta.numericType(N)) {
        .bool => @compileError("zsl.numeric.apery: not defined for " ++ @typeName(N) ++ "."),
        .int => @compileError("zsl.numeric.apery: not defined for " ++ @typeName(N) ++ "."),
        .float => return float.apery(N),
        .dyadic => return dyadic.apery(N),
        .complex => return complex.apery(N), // apery + 0i
        .custom => {
            if (comptime @hasDecl(N, "apery") and @TypeOf(N.apery) == N)
                return N.apery
            else if (comptime meta.hasMethod(N, "apery", fn () N, &.{}))
                return N.apery()
            else
                @compileError("zsl.numeric.apery: " ++ @typeName(N) ++ " must expose a `apery: " ++ @typeName(N) ++ "` declaration, or implement `fn apery() " ++ @typeName(N) ++ "`");
        },
    }
}
