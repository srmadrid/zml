const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const int = @import("../int.zig");

const poly = @import("../poly.zig");

const polutils = @import("utils.zig");

/// Static polynomial type, represented as a contiguous array of
/// `max_degree + 1` coefficients of type `N`.
pub fn Static(max_degree_: comptime_int, N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.poly.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [max_degree + 1]N,

        // Type signatures
        pub const is_poly = true;
        pub const is_static = true;
        pub const max_degree = max_degree_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: poly.Static(max_degree, N) = .{
            .data = undefined,
        };

        pub const zero: poly.Static(max_degree, N) = .{
            .data = @splat(numeric.zero(N)),
        };

        pub const init = empty;

        /// Initializes a new `poly.Static(max_degree, N)`, with all
        /// coefficients set by calling the specified function with the given
        /// arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the polynomial.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `poly.Static(max_degree, N)`: The newly initialized polynomial.
        ///
        /// ## Errors
        pub fn initFn(comptime @"fn": anytype, args: anytype) !poly.Static(max_degree, N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.poly.Static(max_degree, N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var pol: poly.Static(max_degree, N) = .init;

            inline for (0..max_degree + 1) |i| {
                pol.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return pol;
        }

        /// Gets the coefficient of `x^exp`.
        ///
        /// ## Arguments
        /// * `self` (`poly.Static(max_degree, N)`): The polynomial to get the
        ///   coefficient from.
        /// * `exp` (`usize`): The exponent whose coefficient to get.
        ///
        /// ## Returns
        /// `N`: The coefficient of `x^exp`.
        ///
        /// ## Errors
        /// * `poly.Error.PositionOutOfBounds`: If `exp` is out of bounds.
        pub fn get(self: poly.Static(max_degree, N), exp: usize) !N {
            if (exp > max_degree)
                return poly.Error.PositionOutOfBounds;

            return self.getAssumeInBounds(exp);
        }

        /// Gets the coefficient of `x^exp` without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`poly.Static(max_degree, N)`): The polynomial to get the
        ///   coefficient from.
        /// * `exp` (`usize`): The exponent whose coefficient to get.
        ///   Assumed to be within bounds.
        ///
        /// ## Returns
        /// `N`: The coefficient of `x^exp`.
        pub fn getAssumeInBounds(self: poly.Static(max_degree, N), exp: usize) N {
            return self.data[exp];
        }

        /// Sets the coefficient of `x^exp`.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Static(max_degree, N)`): A pointer to the
        ///   polynomial to set the coefficient in.
        /// * `exp` (`usize`): The exponent whose coefficient to set.
        /// * `value` (`N`): The value to set the coefficient to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `poly.Error.PositionOutOfBounds`: If `exp` is out of bounds.
        pub fn set(self: *poly.Static(max_degree, N), exp: usize, value: N) !void {
            if (exp > max_degree)
                return poly.Error.PositionOutOfBounds;

            self.setAssumeInBounds(exp, value);
        }

        /// Sets the coefficient of `x^exp` without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Static(max_degree, N)`): A pointer to the
        ///   polynomial to set the coefficient in.
        /// * `exp` (`usize`): The exponent whose coefficient to set.
        /// * `value` (`N`): The value to set the coefficient to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *poly.Static(max_degree, N), exp: usize, value: N) void {
            self.data[exp] = value;
        }

        /// Returns the degree of the polynomial, i.e., the highest exponent
        /// with a nonzero coefficient.
        ///
        /// ## Arguments
        /// * `self` (`poly.Static(max_degree, N)`): The polynomial to get the
        ///   degree of.
        ///
        /// ## Returns
        /// `void`
        pub fn degree(self: poly.Static(max_degree, N)) usize {
            var i: usize = max_degree + 1;
            while (i > 0) {
                i -= 1;

                if (!numeric.eq(self.data[i], numeric.zero(N)))
                    return i;
            }

            return 0;
        }

        /// Evaluates the polynomial at `x` using Horner's method.
        ///
        /// ## Arguments
        /// * `self` (`poly.Static(max_degree, N)`): The polynomial to evaluate.
        /// * `x` (`N`): The point to evaluate at.
        ///
        /// ## Returns
        /// `N`: `p(x)`.
        pub fn eval(self: poly.Static(max_degree, N), x: N) N {
            var result: meta.Accumulator(N) = numeric.cast(self.data[max_degree]);

            var i: usize = max_degree;
            while (i > 0) {
                i -= 1;

                numeric.addInto(
                    &result,
                    numeric.mul(result, x),
                    self.data[i],
                );
            }

            return numeric.cast(N, result);
        }

        /// Returns the derivative `p'(x)` of a polynomial.
        ///
        /// ## Arguments
        /// * `self` (`poly.Static(max_degree, N)`): The polynomial to get the
        ///   derivative of.
        ///
        /// ## Returns
        /// `poly.Static(max_degree -| 1, N)`: `p'(x)`.
        pub fn derivative(self: poly.Static(max_degree, N)) poly.Static(int.max(1, max_degree -| 1), N) {
            var result: Static(int.max(1, max_degree -| 1), N) = .zero;

            var i: usize = 1;
            while (i <= max_degree) : (i += 1) {
                numeric.mulInto(
                    &result.data[i - 1],
                    self.data[i],
                    i,
                );
            }

            return result;
        }

        /// Returns the antiderivative `c + ∫ p dx` of a polynomial, with the
        /// given constant of integration.
        ///
        /// ## Arguments
        /// * `self` (`poly.Static(max_degree, N)`): The polynomial to
        ///   integrate.
        /// * `constant` (`N`): The constant term of the result.
        ///
        /// ## Returns
        /// `poly.Static(max_degree + 1, N)`: `c + ∫ p dx`.
        pub fn integral(self: poly.Static(max_degree, N), constant: N) Static(max_degree + 1, N) {
            var result: Static(max_degree + 1, N) = .zero;
            result.data[0] = constant;

            var i: usize = 0;
            while (i <= max_degree) : (i += 1) {
                numeric.divInto(
                    &result.data[i + 1],
                    self.data[i],
                    i + 1,
                );
            }

            return result;
        }

        pub fn format(self: poly.Static(max_degree, N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const poly.Static(max_degree, N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .poly = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                poly: *const poly.Static(max_degree, N),

                pub fn format(self: poly.Static(max_degree, N).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const deg = self.poly.degree();

                    try writer.print("zsl.poly.Static({d}, {s}) (degree {d}):\n\n", .{ max_degree, @typeName(N), deg });

                    return polutils.format(self, num_fmt, deg, writer);
                }
            };
        }
    };
}
