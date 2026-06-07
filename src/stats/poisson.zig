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

            return k - 1;
        }
    };
}
