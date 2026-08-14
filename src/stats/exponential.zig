const std = @import("std");

const meta = @import("../meta.zig");

const numeric = @import("../numeric.zig");

const stats = @import("../stats.zig");

const utils = @import("utils.zig");

/// An exponential distribution that yields continuous waiting times of type
/// `N`. The distribution is parameterized by a positive rate `lambda`, it
/// models the interval between independent events occurring continuously at a
/// constant average rate. For complex types, it models independent exponential
/// distributions for the real and imaginary parts.
pub fn Exponential(comptime N: type) type {
    comptime if (!meta.isNumeric(N) or !meta.isNonIntegral(N))
        @compileError("zsl.stats.Exponential: N must be a non-integral numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        lambda: N,

        // Type signatures
        pub const is_distribution = true;

        // Numeric type
        pub const Numeric = N;

        const Self = @This();

        /// Initializes a new exponential distribution.
        ///
        /// ## Arguments
        /// * `lambda` (`N`): The rate parameter (must be positive non-zero).
        pub fn init(lambda: N) Self {
            return .{
                .lambda = lambda,
            };
        }

        /// Samples a random value from the exponential distribution using
        /// Inverse Transform Sampling.
        ///
        /// ## Arguments
        /// * `self` (`stats.Exponential(N)`): The exponential distribution.
        /// * `prng` (`std.Random`): The standard random number generator.
        ///
        /// ## Returns
        /// `N`: A random non-negative waiting time.
        pub fn sample(self: Self, prng: std.Random) N {
            return if (comptime !meta.isComplex(N))
                sampleReal(self.lambda, prng)
            else
                .{
                    .re = sampleReal(self.lambda.re, prng),
                    .im = sampleReal(self.lambda.im, prng),
                };
        }

        fn sampleReal(rate: meta.Real(N), prng: std.Random) meta.Real(N) {
            const u = utils.standardUniform(meta.Real(N), prng);
            const one_minus_u = numeric.sub(1, u);
            return numeric.div(numeric.neg(numeric.ln(one_minus_u)), rate);
        }

        /// Computes the Probability Density Function (PDF) evaluated at `x`.
        /// For complex types, evaluates the joint probability density of the
        /// independent real and imaginary components.
        ///
        /// ## Arguments
        /// * `self` (`stats.Exponential(N)`): The exponential distribution.
        /// * `x` (`N`): The value at which to evaluate the density.
        ///
        /// ## Returns
        /// `meta.Real(N)`: The probability density at `x`.
        pub fn pdf(self: Self, x: N) meta.Real(N) {
            return if (comptime !meta.isComplex(N))
                pdfReal(x, self.lambda)
            else
                numeric.mul(
                    pdfReal(x.re, self.lambda.re),
                    pdfReal(x.im, self.lambda.im),
                );
        }

        fn pdfReal(val: meta.Real(N), rate: meta.Real(N)) meta.Real(N) {
            if (numeric.lt(val, 0))
                return numeric.cast(meta.Real(N), 0);

            const exp_term = numeric.exp(numeric.neg(numeric.mul(rate, val)));
            return numeric.mul(rate, exp_term);
        }

        /// Computes the natural logarithm of the Probability Density Function
        /// (Log-PDF) evaluated at `x`. Recommended over `ln(pdf(x))` to prevent
        /// underflow for large waiting times.
        ///
        /// ## Arguments
        /// * `self` (`stats.Exponential(N)`): The exponential distribution.
        /// * `x` (`N`): The value at which to evaluate the log-density.
        ///
        /// ## Returns
        /// `meta.Real(N)`: The log-probability density at `x`.
        pub fn lpdf(self: Self, x: N) meta.Real(N) {
            return if (comptime !meta.isComplex(N))
                lpdfReal(x, self.lambda)
            else
                numeric.add(
                    lpdfReal(x.re, self.lambda.re),
                    lpdfReal(x.im, self.lambda.im),
                );
        }

        fn lpdfReal(val: meta.Real(N), rate: meta.Real(N)) meta.Real(N) {
            if (numeric.lt(val, 0))
                return numeric.neg(numeric.inf(meta.Real(N)));

            // lpdf(x) = ln(lambda) - lambda * x
            const ln_rate = numeric.ln(rate);
            const rate_x = numeric.mul(rate, val);
            return numeric.sub(ln_rate, rate_x);
        }

        /// Computes the Cumulative Distribution Function (CDF) evaluated at
        /// `x`. Represents the probability that an event occurs within waiting
        /// time `x`.
        ///
        /// ## Arguments
        /// * `self` (`stats.Exponential(N)`): The exponential distribution.
        /// * `x` (`N`): The upper bound of the waiting time.
        ///
        /// ## Returns
        /// `meta.Real(N)`: The cumulative probability in the range [0, 1].
        pub fn cdf(self: Self, x: N) meta.Real(N) {
            return if (comptime !meta.isComplex(N))
                cdfReal(x, self.lambda)
            else
                numeric.mul(
                    cdfReal(x.re, self.lambda.re),
                    cdfReal(x.im, self.lambda.im),
                );
        }

        fn cdfReal(val: meta.Real(N), rate: meta.Real(N)) meta.Real(N) {
            if (numeric.le(val, 0))
                return numeric.cast(meta.Real(N), 0);

            // cdf(x) = 1 - e^(-lambda * x)
            const exp_term = numeric.exp(numeric.neg(numeric.mul(rate, val)));
            return numeric.sub(1, exp_term);
        }

        /// Computes the Inverse Cumulative Distribution Function (iCDF), or
        /// quantile function. Maps a probability threshold `p` to the
        /// corresponding waiting time `x`.
        ///
        /// ## Arguments
        /// * `self` (`stats.Exponential(N)`): The exponential distribution.
        /// * `p` (`N`): The probability threshold. Must be in [0, 1) for real
        ///   types, or have components in [0, 1) for complex types.
        ///
        /// ## Returns
        /// `N`: The waiting time `x` such that `cdf(x) == p`.
        pub fn icdf(self: Self, p: N) N {
            return if (comptime !meta.isComplex(N))
                icdfReal(p, self.lambda)
            else
                .{
                    .re = icdfReal(p.re, self.lambda.re),
                    .im = icdfReal(p.im, self.lambda.im),
                };
        }

        fn icdfReal(prob: meta.Real(N), rate: meta.Real(N)) meta.Real(N) {
            if (numeric.le(prob, 0))
                return numeric.cast(meta.Real(N), 0);

            if (numeric.ge(prob, 1))
                return numeric.inf(meta.Real(N));

            // icdf(p) = -ln(1 - p) / lambda
            const one_minus_p = numeric.sub(1, prob);
            return numeric.div(numeric.neg(numeric.ln(one_minus_p)), rate);
        }
    };
}
