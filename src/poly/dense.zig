const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const int = @import("../int.zig");

const poly = @import("../poly.zig");

const polyutils = @import("utils.zig");

/// Dense polynomial type, represented as a contiguous array of `max_degree + 1`
/// coefficients of type `N`.
pub fn Dense(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.poly.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        max_degree: usize,
        flags: poly.Flags,

        // Type signatures
        pub const is_poly = true;
        pub const is_dense = true;

        // Numeric type
        pub const Numeric = N;

        pub const empty = Dense(N){
            .data = &.{},
            .max_degree = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `poly.Dense(N)` with `max_degree` coefficients.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `max_degree` (`usize`): The maximum number of coefficients.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: The newly initialized polynomial.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn init(allocator: std.mem.Allocator, max_degree: usize) !poly.Dense(N) {
            return .{
                .data = (try allocator.alloc(N, max_degree + 1)).ptr,
                .max_degree = max_degree,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new `poly.Dense(N)` with the given buffer.
        ///
        /// ## Arguments
        /// * `buffer` (`[]N`): The buffer.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: The newly initialized polynomial.
        ///
        /// ## Errors
        /// * `poly.Error.ZeroLength`: If if the length of the buffer is zero.
        pub fn initBuffer(buffer: []N) !poly.Dense(N) {
            if (buffer.len == 0)
                return poly.Error.ZeroLength;

            return .{
                .data = buffer.ptr,
                .max_degree = buffer.len - 1,
                .flags = .{ .owns_data = false },
            };
        }

        /// Initializes a new `poly.Dense(N)` with `max_degree` coefficients,
        /// all set to the specialized value.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `max_degree` (`usize`): The number of coefficients.
        /// * `value` (`N`): The value to fill the polynomial with.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: The newly initialized polynomial.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.ZeroMaxDegree`: If `max_degree` is zero.
        pub fn initValue(allocator: std.mem.Allocator, max_degree: usize, value: N) !poly.Dense(N) {
            const pol: poly.Dense(N) = try .init(allocator, max_degree);
            @memset(pol.data[0 .. max_degree + 1], value);

            return pol;
        }

        /// Initializes a new `poly.Dense(N)` with `max_degree` coefficients,
        /// with all coefficients set by calling the specified function with the
        /// given arguments.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `max_degree` (`usize`): The number of coefficients.
        /// * `@"fn"` (`anytype`): The function to call to fill the polynomial.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `poly.Static(max_degree, N)`: The newly initialized polynomial.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn initFn(allocator: std.mem.Allocator, max_degree: usize, comptime @"fn": anytype, args: anytype) !poly.Dense(N) {
            var pol: poly.Dense(N) = try .init(allocator, max_degree);
            errdefer pol.deinit(allocator);

            var i: usize = 0;
            while (i <= max_degree) : (i += 1) {
                pol.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(@TypeOf(.{i} ++ args)))) == .error_union)
                    try @call(.auto, @"fn", .{i} ++ args)
                else
                    @call(.auto, @"fn", .{i} ++ args);
            }

            return pol;
        }

        /// Initializes the monic polynomial
        /// `(x - r[0])(x - r[1])⋯(x - r[n - 1])` from its roots.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `roots` (`[]const N`): The roots of the desired polynomial.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: The polynomial with the given roots.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation
        ///   fails.
        /// * `poly.Error.ZeroLength`: If if the length of the roots buffer is
        ///   zero.
        pub fn initRoots(allocator: std.mem.Allocator, roots: []const N) !poly.Dense(N) {
            if (roots.len == 0)
                return poly.Error.ZeroLength;

            var pol: poly.Dense(N) = try .initValue(allocator, roots.len, numeric.zero(N));

            pol.data[0] = numeric.one(N);

            for (roots, 0..) |r, k| {
                var j: usize = k + 1;
                while (j > 0) {
                    numeric.subInto(
                        &pol.data[j],
                        pol.data[j - 1],
                        numeric.mul(pol.data[j], r),
                    );

                    j -= 1;
                }

                numeric.subInto(
                    &pol.data[0],
                    numeric.zero(N),
                    numeric.mul(pol.data[0], r),
                );
            }

            return pol;
        }

        /// Deinitializes the polynomial, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Dense(N)`): A pointer to the polynomial to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *poly.Dense(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.max_degree + 1]);
            }

            self.* = undefined;
        }

        /// Gets the coefficient of `x^exp`.
        ///
        /// ## Arguments
        /// * `self` (`poly.Dense(N)`): The polynomial to get the coefficient
        ///   from.
        /// * `exp` (`usize`): The exponent whose coefficient to get.
        ///
        /// ## Returns
        /// `N`: The coefficient of `x^exp`.
        ///
        /// ## Errors
        /// * `poly.Error.PositionOutOfBounds`: If `exp` is out of bounds.
        pub fn get(self: poly.Dense(N), exp: usize) !N {
            if (exp > self.max_degree)
                return poly.Error.PositionOutOfBounds;

            return self.getAssumeInBounds(exp);
        }

        /// Gets the coefficient of `x^exp` without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`poly.Dense(N)`): The polynomial to get the coefficient
        ///   from.
        /// * `exp` (`usize`): The exponent whose coefficient to get.
        ///   Assumed to be within bounds.
        ///
        /// ## Returns
        /// `N`: The coefficient of `x^exp`.
        pub fn getAssumeInBounds(self: poly.Dense(N), exp: usize) N {
            return self.data[exp];
        }

        /// Sets the coefficient of `x^exp`.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Dense(N)`): A pointer to the polynomial to set the
        ///   coefficient in.
        /// * `exp` (`usize`): The exponent whose coefficient to set.
        /// * `value` (`N`): The value to set the coefficient to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `poly.Error.PositionOutOfBounds`: If `exp` is out of bounds.
        pub fn set(self: *poly.Dense(N), exp: usize, value: N) !void {
            if (exp > self.max_degree)
                return poly.Error.PositionOutOfBounds;

            self.setAssumeInBounds(exp, value);
        }

        /// Sets the coefficient of `x^exp` without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Dense(N)`): A pointer to the polynomial to set the
        ///   coefficient in.
        /// * `exp` (`usize`): The exponent whose coefficient to set.
        /// * `value` (`N`): The value to set the coefficient to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *poly.Dense(N), exp: usize, value: N) void {
            self.data[exp] = value;
        }

        /// Reserves space for at least `max_degree` coefficients, zeroing out
        /// any extra allocated coefficients.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Dense(N)`): A pointer to the polynomial to reserve
        ///   space for.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `new_max_degree` (`usize`): The new capacity for coefficients.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.ZeroMaxDegree`: If `new_max_degree` is zero.
        pub fn reserve(self: *poly.Dense(N), allocator: std.mem.Allocator, new_max_degree: usize) !void {
            if (!self.flags.owns_data or new_max_degree <= self.max_degree)
                return;

            if (new_max_degree == 0)
                return poly.Error.ZeroMaxDegree;

            const new_data = try allocator.realloc(self.data[0 .. self.max_degree + 1], new_max_degree + 1);

            if (new_max_degree > self.max_degree)
                @memset(new_data[self.max_degree..new_max_degree], numeric.zero(N));

            self.data = new_data.ptr;
            self.max_degree = new_max_degree;
        }

        /// Creates a copy of the polynomial.
        ///
        /// ## Arguments
        /// * `self` (`poly.Dense(N)`): The polynomial to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: The copied polynomial.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: poly.Dense(N), allocator: std.mem.Allocator) !poly.Dense(N) {
            var pol: poly.Dense(N) = try .init(allocator, self.max_degree);

            for (0..pol.max_degree + 1) |i| {
                pol.data[i] = self.data[i];
            }

            return pol;
        }

        /// Returns the degree of the polynomial, i.e., the highest exponent
        /// with a nonzero coefficient.
        ///
        /// ## Arguments
        /// * `self` (`poly.Dense(N)`): The polynomial to get the degree of.
        ///
        /// ## Returns
        /// `void`
        pub fn degree(self: poly.Dense(N)) usize {
            var i: usize = self.max_degree + 1;
            while (i > 0) {
                i -= 1;

                if (!numeric.eq(self.data[i], numeric.zero(N)))
                    return i;
            }

            return 0;
        }

        /// Returns a copy with trailing zero coefficients removed, so that
        /// `result.max_degree == result.degree()`.
        ///
        /// ## Arguments
        /// * `self` (`poly.Dense(N)`): The polynomial to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: The copied polynomial.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn trimmed(self: poly.Dense(N), allocator: std.mem.Allocator) !poly.Dense(N) {
            const deg = self.degree();

            var pol: poly.Dense(N) = try .init(allocator, deg);

            for (0..deg + 1) |i| {
                pol.data[i] = self.data[i];
            }

            return pol;
        }

        /// Evaluates the polynomial at `x` using Horner's method.
        ///
        /// ## Arguments
        /// * `self` (`poly.Dense(N)`): The polynomial to evaluate.
        /// * `x` (`N`): The point to evaluate at.
        ///
        /// ## Returns
        /// `N`: `p(x)`.
        pub fn eval(self: poly.Dense(N), x: N) N {
            var result: meta.Accumulator(N) = numeric.cast(self.data[self.max_degree]);

            var i: usize = self.max_degree;
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
        /// * `self` (`poly.Dense(N)`): The polynomial to get the derivative of.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: `p'(x)`.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn derivative(self: poly.Dense(N), allocator: std.mem.Allocator) !poly.Dense(N) {
            if (self.max_degree == 0)
                return .initValue(allocator, 1, numeric.zero(N));

            var result: Dense(N) = try .init(allocator, int.max(1, self.max_degree -| 1));

            var i: usize = 1;
            while (i <= self.max_degree) : (i += 1) {
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
        /// * `self` (`poly.Dense(N)`): The polynomial to integrate.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `constant` (`N`): The constant term of the result.
        ///
        /// ## Returns
        /// `poly.Dense(N)`: `c + ∫ p dx`.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn integral(self: poly.Dense(N), allocator: std.mem.Allocator, constant: N) !Dense(N) {
            var result: Dense(N) = try .init(allocator, self.max_degree + 1);
            result.data[0] = constant;

            var i: usize = 0;
            while (i <= self.max_degree) : (i += 1) {
                numeric.divInto(
                    &result.data[i + 1],
                    self.data[i],
                    i + 1,
                );
            }

            return result;
        }

        pub fn format(self: poly.Dense(N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const poly.Dense(N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .poly = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                poly: *const poly.Dense(N),

                pub fn format(self: Dense(N).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const deg = self.poly.degree();

                    try writer.print("zsl.poly.Dense({s}) (degree {d}):\n\n", .{ @typeName(N), deg });

                    return polyutils.format(self, num_fmt, deg, writer);
                }
            };
        }
    };
}
