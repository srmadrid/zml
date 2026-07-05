const std = @import("std");

const meta = @import("../meta.zig");

const numeric = @import("../numeric.zig");

const stats = @import("../stats.zig");

const utils = @import("utils.zig");

/// A Poisson distribution that yields discrete event counts of type `Int`. The
/// distribution is parameterized by a continuous rate `lambda` of type `Real`.
pub fn Poisson(comptime Int: type, comptime Real: type) type {
    comptime if (!meta.isNumeric(Int) or meta.isNonIntegral(Int) or !meta.isNumeric(Real) or meta.isIntegral(Real) or !meta.isReal(Real))
        @compileError("zsl.stats.Poisson: Int must be an integral numeric type, and Real must be a real non-integral numeric type, got \n\tInt = " ++ @typeName(Int) ++ "\n\tReal = " ++ @typeName(Real) ++ "\n");

    return struct {
        lambda: Real,
        enl: Real, // Cached e^-λ

        // Type signatures
        pub const is_distribution = true;

        // Numeric type
        pub const Numeric = Int;

        /// Initializes a new Poisson distribution.
        ///
        /// ## Arguments
        /// * `lambda` (`Real`): The expected rate of occurrences (must be
        ///   positive non-zero).
        pub fn init(lambda: Real) Poisson(Int, Real) {
            return .{
                .lambda = lambda,
                .enl = numeric.exp(numeric.neg(lambda)),
            };
        }

        /// Samples a random integer from the Poisson distribution using Knuth's
        /// multiplicative method.
        ///
        /// ## Arguments
        /// * `self` (`stats.Poisson(Int, Real)`): The Poisson distribution to
        ///   sample from.
        /// * `prng` (`std.Random`): The standard random number generator.
        ///
        /// ## Returns
        /// `Int`: A random number of events generated based on `lambda`.
        pub fn sample(self: Poisson(Int, Real), prng: std.Random) Int {
            var k = numeric.zero(Int);
            var p = numeric.one(Real);

            while (numeric.gt(p, self.enl)) {
                numeric.addInto(&k, k, numeric.one(Int));
                const u = utils.standardUniform(Real, prng);
                numeric.mulInto(&p, p, u);
            }

            return numeric.sub(k, numeric.one(Int));
        }

        /// Computes the Probability Mass Function (PMF) evaluated at `k`.
        /// Represents the exact probability of observing `k` events in a given
        /// interval with expected rate `lambda`.
        ///
        /// ## Arguments
        /// * `self` (`stats.Poisson(Int, Real)`): The Poisson distribution.
        /// * `k` (`Int`): The number of occurrences to evaluate.
        ///
        /// ## Returns
        /// `Real`: The probability mass at `k` in the range [0, 1].
        pub fn pmf(self: Poisson(Int, Real), k: Int) Real {
            if (numeric.lt(k, numeric.zero(Int)))
                return numeric.zero(Real);

            return numeric.exp(self.lpmf(k));
        }

        /// Computes the natural logarithm of the Probability Mass Function
        /// (Log-PMF) evaluated at `k`. Recommended over `ln(pmf(k))` for large
        /// values of `k` or `lambda` to prevent floating-point underflow.
        ///
        /// ## Arguments
        /// * `self` (`stats.Poisson(Int, Real)`): The Poisson distribution.
        /// * `k` (`Int`): The number of occurrences to evaluate.
        ///
        /// ## Returns
        /// `Real`: The natural logarithm of the probability mass at `k`.
        pub fn lpmf(self: Poisson(Int, Real), k: Int) Real {
            if (numeric.lt(k, numeric.zero(Int)))
                return numeric.neg(numeric.inf(Real));

            // lpmf(k) = k * ln(lambda) - lambda - ln(k!), with k! = lgamma(k + 1)
            const k_ln_lambda = numeric.mul(numeric.cast(Real, k), numeric.ln(self.lambda));
            const log_factorial = numeric.lgamma(numeric.add(numeric.cast(Real, k), numeric.one(Real)));

            const diff = numeric.sub(k_ln_lambda, self.lambda);
            return numeric.sub(diff, log_factorial);
        }

        /// Computes the Cumulative Distribution Function (CDF) evaluated at
        /// `k`. Represents the probability of observing `k` or fewer events.
        ///
        /// Uses an iterative recurrence relation to sum probability masses
        /// efficiently without redundant power or factorial calculations.
        ///
        /// ## Arguments
        /// * `self` (`stats.Poisson(Int, Real)`): The Poisson distribution.
        /// * `k` (`Int`): The upper bound of occurrences to include.
        ///
        /// ## Returns
        /// `Real`: The cumulative probability in the range [0, 1].
        pub fn cdf(self: Poisson(Int, Real), k: Int) Real {
            if (numeric.lt(k, numeric.zero(Int)))
                return numeric.zero(Real);

            var sum = self.enl; // P(X = 0) = e^-lambda
            var term = self.enl;
            var i: Int = numeric.one(Int);

            while (numeric.le(i, k)) : (numeric.addInto(&i, numeric.one(Int))) {
                // P(X = i) = P(X = i - 1) * (lambda / i)
                term = numeric.mul(term, numeric.div(self.lambda, numeric.cast(Real, i)));
                sum = numeric.add(sum, term);
            }

            if (numeric.gt(sum, numeric.one(Real)))
                return numeric.one(Real);

            return sum;
        }

        /// Computes the Inverse Cumulative Distribution Function (iCDF), or
        /// quantile function. For a discrete Poisson distribution, this returns
        /// the smallest integer `k` such that `cdf(k) >= p`.
        ///
        /// ## Arguments
        /// * `self` (`stats.Poisson(Int, Real)`): The Poisson distribution.
        /// * `p` (`Real`): The probability threshold in the range [0, 1].
        ///
        /// ## Returns
        /// `Int`: The smallest number of occurrences whose cumulative
        /// probability meets or exceeds `p`.
        pub fn icdf(self: Poisson(Int, Real), p: Real) Int {
            if (numeric.le(p, numeric.zero(Real)))
                return 0;

            if (numeric.ge(p, numeric.one(Real)))
                return numeric.maxVal(Int);

            var k = numeric.zero(Int);
            var term = self.enl; // P(X = 0)
            var sum = term;

            while (numeric.lt(sum, p)) {
                numeric.addInto(&k, numeric.one(Int));
                const k_real = numeric.cast(Real, k);
                term = numeric.mul(term, numeric.div(self.lambda, k_real));
                sum = numeric.add(sum, term);
            }

            return k;
        }
    };
}
