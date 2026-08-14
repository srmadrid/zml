const std = @import("std");

const meta = @import("../meta.zig");

const numeric = @import("../numeric.zig");

const stats = @import("../stats.zig");

const utils = @import("utils.zig");

/// A normal (Gaussian) distribution that yields values of type `N`.
/// For complex types, the distribution generates independent normal variables
/// for the real and imaginary parts, scaled by `sigma.re` and `sigma.im`
/// respectively. If `sigma.re == sigma.im`, this results in a circularly
/// symmetric complex normal distribution.
pub fn Normal(comptime N: type) type {
    comptime if (!meta.isNumeric(N) or !meta.isNonIntegral(N))
        @compileError("zsl.stats.Normal: N must be a non-integral numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        mu: N,
        sigma: N,

        // Type signatures
        pub const is_distribution = true;

        // Numeric type
        pub const Numeric = N;

        const Self = @This();

        /// Initializes a new normal distribution.
        ///
        /// ## Arguments
        /// * `mu` (`N`): The mean (expected value) of the distribution.
        /// * `sigma` (`N`): The standard deviation of the distribution.
        pub fn init(mu: N, sigma: N) Self {
            return .{
                .mu = mu,
                .sigma = sigma,
            };
        }

        /// Samples a random value from the normal distribution using the
        /// Marsaglia polar method.
        ///
        /// ## Arguments
        /// * `self` (`stats.Normal(N)`): The normal distribution to sample
        ///   from.
        /// * `prng` (`std.Random`): The standard random number generator used
        ///   to produce the random bits.
        ///
        /// ## Returns
        /// `N`: A random value normally distributed with mean `mu` and standard
        /// deviation `sigma`.
        pub fn sample(self: Self, prng: std.Random) N {
            var u: meta.Real(N) = undefined;
            var v: meta.Real(N) = undefined;
            var s: meta.Real(N) = undefined;

            while (true) {
                u = numeric.sub(numeric.mul(2, utils.standardUniform(meta.Real(N), prng)), 1);
                v = numeric.sub(numeric.mul(2, utils.standardUniform(meta.Real(N), prng)), 1);
                s = numeric.add(numeric.mul(u, u), numeric.mul(v, v));

                if (numeric.gt(s, 0) and numeric.lt(s, 1))
                    break;
            }

            const temp = numeric.sqrt(numeric.mul(-2, numeric.div(numeric.ln(s), s)));

            return if (comptime !meta.isComplex(N))
                numeric.add(self.mu, numeric.mul(self.sigma, numeric.mul(u, temp)))
            else
                .{
                    .re = numeric.add(self.mu.re, numeric.mul(self.sigma.re, numeric.mul(u, temp))),
                    .im = numeric.add(self.mu.im, numeric.mul(self.sigma.im, numeric.mul(v, temp))),
                };
        }

        /// Computes the Probability Density Function (PDF) evaluated at `x`.
        /// For complex types, this evaluates the joint probability density of
        /// the independent real and imaginary components.
        ///
        /// ## Arguments
        /// * `self` (`stats.Normal(N)`): The normal distribution.
        /// * `x` (`N`): The value at which to evaluate the density.
        ///
        /// ## Returns
        /// `meta.Real(N)`: The probability density at `x`.
        pub fn pdf(self: Self, x: N) meta.Real(N) {
            return if (comptime !meta.isComplex(N))
                pdfReal(x, self.mu, self.sigma)
            else
                numeric.mul(
                    pdfReal(x.re, self.mu.re, self.sigma.re),
                    pdfReal(x.im, self.mu.im, self.sigma.im),
                );
        }

        fn pdfReal(val: meta.Real(N), mean: meta.Real(N), std_dev: meta.Real(N)) meta.Real(N) {
            const z = numeric.div(numeric.sub(val, mean), std_dev);
            const z_sq = numeric.mul(z, z);
            const half_z_sq = numeric.div(z_sq, 2);
            const exp_term = numeric.exp(numeric.neg(half_z_sq));

            const sqrt_two_pi = numeric.sqrt(numeric.mul(2, numeric.pi(meta.Real(N))));
            const norm_factor = numeric.mul(std_dev, sqrt_two_pi);

            return numeric.div(exp_term, norm_factor);
        }

        /// Computes the natural logarithm of the Probability Density Function
        /// (Log-PDF) evaluated at `x`. This is numerically much more stable
        /// than taking the log of `pdf(x)` for extreme values.
        ///
        /// ## Arguments
        /// * `self` (`stats.Normal(N)`): The normal distribution.
        /// * `x` (`N`): The value at which to evaluate the log-density.
        ///
        /// ## Returns
        /// `meta.Real(N)`: The log-probability density at `x`.
        pub fn lpdf(self: Self, x: N) meta.Real(N) {
            return if (comptime !meta.isComplex(N))
                lpdfReal(x, self.mu, self.sigma)
            else
                numeric.add(
                    lpdfReal(x.re, self.mu.re, self.sigma.re),
                    lpdfReal(x.im, self.mu.im, self.sigma.im),
                );
        }

        fn lpdfReal(val: meta.Real(N), mean: meta.Real(N), std_dev: meta.Real(N)) meta.Real(N) {
            const z = numeric.div(numeric.sub(val, mean), std_dev);
            const z_sq = numeric.mul(z, z);
            const half_z_sq = numeric.div(z_sq, 2);

            const log_sigma = numeric.ln(std_dev);
            const log_two_pi = numeric.ln(numeric.mul(2, numeric.pi(meta.Real(N))));
            const half_log_two_pi = numeric.div(log_two_pi, 2);

            const norm_term = numeric.add(log_sigma, half_log_two_pi);
            return numeric.neg(numeric.add(norm_term, half_z_sq));
        }

        /// Computes the Cumulative Distribution Function (CDF) evaluated at
        /// `x`. Represents the probability that a random variable from this
        /// distribution is less than or equal to `x`. For complex types, this
        /// returns the joint cumulative probability of both real and imaginary
        /// components.
        ///
        /// ## Arguments
        /// * `self` (`stats.Normal(N)`): The normal distribution.
        /// * `x` (`N`): The upper bound to evaluate the cumulative probability.
        ///
        /// ## Returns
        /// `meta.Real(N)`: The cumulative probability in the range [0, 1].
        pub fn cdf(self: Self, x: N) meta.Real(N) {
            return if (comptime !meta.isComplex(N))
                cdfReal(x, self.mu, self.sigma)
            else
                numeric.mul(
                    cdfReal(x.re, self.mu.re, self.sigma.re),
                    cdfReal(x.im, self.mu.im, self.sigma.im),
                );
        }

        fn cdfReal(val: meta.Real(N), mean: meta.Real(N), std_dev: meta.Real(N)) meta.Real(N) {
            const diff = numeric.sub(val, mean);
            const denom = numeric.mul(std_dev, numeric.sqrt(numeric.cast(meta.Real(N), 2)));
            const erf_val = numeric.erf(numeric.div(diff, denom));

            const one_plus_erf = numeric.add(1, erf_val);
            return numeric.div(one_plus_erf, 2);
        }

        /// Computes the Inverse Cumulative Distribution Function (iCDF), also
        /// known as the quantile function or percent point function. Maps a
        /// probability `p` back to the corresponding value in the distribution
        /// domain.
        ///
        /// ## Arguments
        /// * `self` (`stats.Normal(N)`): The normal distribution.
        /// * `p` (`N`): The probability threshold. For real types, p must be in
        ///   (0, 1). For complex types, both p.re and p.im must be in (0, 1).
        ///
        /// ## Returns
        /// `N`: The value `x` such that `cdf(x) == p`.
        pub fn icdf(self: Self, p: N) N {
            return if (comptime !meta.isComplex(N))
                icdfReal(p, self.mu, self.sigma)
            else
                .{
                    .re = icdfReal(p.re, self.mu.re, self.sigma.re),
                    .im = icdfReal(p.im, self.mu.im, self.sigma.im),
                };
        }

        fn icdfReal(prob: meta.Real(N), mean: meta.Real(N), std_dev: meta.Real(N)) meta.Real(N) {
            const two_p_minus_one = numeric.sub(numeric.mul(2, prob), 1);
            const erf_inv_val = numeric.erfinv(two_p_minus_one);

            const sqrt_two = numeric.sqrt(numeric.cast(meta.Real(N), 2));
            const scaled = numeric.mul(std_dev, numeric.mul(sqrt_two, erf_inv_val));
            return numeric.add(mean, scaled);
        }
    };
}
