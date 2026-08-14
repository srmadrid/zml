const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const poly = @import("../poly.zig");

const polyutils = @import("utils.zig");

/// Sparse polynomial type, represented as a contiguous array of non-zero
/// coefficients of type `N` along with their corresponding exponents, sorted
/// in ascending order and without duplicates.
pub fn Sparse(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.poly.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        exps: [*]usize,
        nnz: usize,
        _dlen: usize, // allocated length of data
        _elen: usize, // allocated length of exps
        flags: poly.Flags,

        // Type signatures
        pub const is_poly = true;
        pub const is_sparse = true;

        // Numeric type
        pub const Numeric = N;

        pub const empty: poly.Sparse(N) = .{
            .data = &.{},
            .exps = &.{},
            .nnz = 0,
            ._dlen = 0,
            ._elen = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new sparse polynomial with an initial capacity for
        /// non-zero coefficients.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `nnz` (`usize`): The initial capacity for non-zero coefficients.
        ///
        /// ## Returns
        /// `vector.Sparse(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn init(allocator: std.mem.Allocator, nnz: usize) !poly.Sparse(N) {
            const data: []N = try allocator.alloc(N, nnz);
            errdefer allocator.free(data);

            return .{
                .data = data.ptr,
                .exps = (try allocator.alloc(usize, nnz)).ptr,
                .nnz = 0,
                ._dlen = nnz,
                ._elen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        /// Deinitializes the polynomial, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Sparse(N)`): A pointer to the polynomial to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *poly.Sparse(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.exps[0..self._elen]);
            }

            self.* = undefined;
        }

        /// Reserves space for at least `new_nnz` non-zero coefficients. If
        /// `self` does not own its data or if `new_nnz` is less than or equal
        /// to the current capacity, this function is a no-op.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Sparse(N)`): A pointer to the polynomial to reserve
        ///   space for.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `new_nnz` (`usize`): The new capacity for non-zero coefficients.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn reserve(self: *poly.Sparse(N), allocator: std.mem.Allocator, new_nnz: usize) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._elen)
                return;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._elen) {
                self.exps = (try allocator.realloc(self.exps[0..self._elen], new_nnz)).ptr;
                self._elen = new_nnz;
            }
        }

        fn find(self: poly.Sparse(N), exp: usize) ?usize {
            var lo: usize = 0;
            var hi: usize = self.nnz;

            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                const mid_exp = self.exps[mid];

                if (mid_exp == exp)
                    return mid;

                if (mid_exp < exp)
                    lo = mid + 1
                else
                    hi = mid;
            }

            return null;
        }

        /// Gets the coefficient of `x^exp`.
        ///
        /// ## Arguments
        /// * `self` (`poly.Sparse(N)`): The polynomial to get the coefficient
        ///   from.
        /// * `exp` (`usize`): The exponent whose coefficient to get.
        ///
        /// ## Returns
        /// `N`: The coefficient of `x^exp`.
        pub fn get(self: poly.Sparse(N), exp: usize) N {
            var lo: usize = 0;
            var hi: usize = self.nnz;

            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                const mid_exp = self.exps[mid];

                if (mid_exp == exp)
                    return self.data[mid];

                if (mid_exp < exp)
                    lo = mid + 1
                else
                    hi = mid;
            }

            return numeric.cast(N, 0);
        }

        /// Sets the element of `x^exp`, inserting it if it does not already
        /// exist and shifting elements as necessary to maintain exponent order.
        ///
        /// ## Arguments
        /// * `self` (`*poly.Sparse(N)`): A pointer to the polynomial to set the
        ///   coefficient in.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `exp` (`usize`): The exponent whose coefficient to set.
        /// * `value` (`N`): The value to set the coefficient to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails
        ///   when inserting a new element.
        /// * `poly.Error.DataNotOwned`: If the polynomial does not own its data
        ///   and a resize is required.
        pub fn set(self: *poly.Sparse(N), allocator: std.mem.Allocator, exp: usize, value: N) !void {
            var lo: usize = 0;
            var hi: usize = self.nnz;

            while (lo < hi) {
                const mid = lo + (hi - lo) / 2;
                const mid_exp = self.exps[mid];

                if (mid_exp == exp) {
                    self.data[mid] = value;
                    return;
                }

                if (mid_exp < exp)
                    lo = mid + 1
                else
                    hi = mid;
            }

            if (self.nnz == self._dlen or self.nnz == self._elen) {
                if (!self.flags.owns_data)
                    return poly.Error.DataNotOwned;

                // Need more space
                var new_nnz = self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to the right to make space for the new element
            const i: usize = lo;
            var j: usize = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.exps[j] = self.exps[j - 1];
            }

            self.data[i] = value;
            self.exps[i] = exp;
            self.nnz += 1;
        }

        /// Returns the degree of the polynomial, i.e., the highest exponent
        /// with a nonzero coefficient.
        ///
        /// ## Arguments
        /// * `self` (`poly.Sparse(N)`): The polynomial to get the degree of.
        ///
        /// ## Returns
        /// `void`
        pub fn degree(self: poly.Sparse(N)) usize {
            return if (self.nnz == 0) 0 else self.exps[self.nnz - 1];
        }

        /// Evaluates the polynomial at `x`.
        ///
        /// ## Arguments
        /// * `self` (`poly.Sparse(N)`): The polynomial to evaluate.
        /// * `x` (`N`): The point to evaluate at.
        ///
        /// ## Returns
        /// `N`: `p(x)`.
        pub fn eval(self: poly.Sparse(N), x: N) N {
            var result: N = numeric.cast(N, 0);

            for (0..self.nnz) |i|
                numeric.addInto(
                    &result,
                    result,
                    numeric.mul(
                        self.data[i],
                        numeric.pow(x, self.exps[i]),
                    ),
                );

            return result;
        }

        /// Returns the derivative `p'(x)` of a polynomial.
        ///
        /// ## Arguments
        /// * `self` (`poly.Sparse(N)`): The polynomial to get the derivative of.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `poly.Sparse(N)`: `p'(x)`.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn derivative(self: poly.Sparse(N), allocator: std.mem.Allocator) !poly.Sparse(N) {
            const skip_constant: usize = if (self.nnz > 0 and self.exps[0] == 0) 1 else 0;
            var result: Sparse(N) = try .init(allocator, self.nnz - skip_constant);

            for (skip_constant..self.nnz, 0..) |i, j| {
                numeric.mulInto(
                    &result.data[j],
                    self.data[i],
                    self.exps[i],
                );

                result.exps[j] = self.exps[i] - 1;
            }

            result.nnz = self.nnz - skip_constant;

            return result;
        }

        /// Returns the antiderivative `c + ∫ p dx` of a polynomial, with the
        /// given constant of integration.
        ///
        /// ## Arguments
        /// * `self` (`poly.Sparse(N)`): The polynomial to integrate.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `constant` (`N`): The constant term of the result.
        ///
        /// ## Returns
        /// `poly.Sparse(N)`: `c + ∫ p dx`.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn integral(self: poly.Sparse(N), allocator: std.mem.Allocator, constant: N) !poly.Sparse(N) {
            const has_constant: usize = if (numeric.eq(constant, numeric.cast(N))) 0 else 1;
            var result: Sparse(N) = try .init(allocator, self.nnz + has_constant);

            if (has_constant == 1) {
                result.data[0] = constant;
                result.exps[0] = 0;
            }

            for (0..self.nnz) |i| {
                numeric.divInto(
                    &result.data[i + has_constant],
                    self.data[i],
                    self.exps[i] + 1,
                );

                result.exps[i + has_constant] = self.exps[i] + 1;
            }

            result.nnz = self.nnz + has_constant;

            return result;
        }

        pub fn format(self: poly.Sparse(N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const poly.Sparse(N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .poly = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                poly: *const poly.Sparse(N),

                pub fn format(self: Sparse(N).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const gap_size: usize = 2;
                    const gap_str = " " ** gap_size ++ "+" ++ " " ** gap_size;

                    try writer.print("zsl.poly.Sparse({s}) (degree {d}):\n\n", .{ @typeName(N), self.poly.degree() });

                    try writer.writeAll("          ");

                    var printed_any_terms = false;

                    var i: usize = self.poly.nnz;
                    while (i > 0) {
                        i -= 1;

                        const val = self.poly.data[i];
                        if (numeric.eq(val, 0))
                            continue;

                        if (printed_any_terms)
                            try writer.writeAll(gap_str);

                        try writer.print(num_fmt, .{val});

                        const exp = self.poly.exps[i];
                        if (exp == 1) {
                            try writer.writeAll(" x");
                        } else if (exp > 1) {
                            try writer.writeAll(" x");
                            try polyutils.writeSuperscript(writer, exp);
                        }

                        printed_any_terms = true;
                    }

                    if (!printed_any_terms)
                        try writer.print(num_fmt, .{numeric.cast(N, 0)});
                }
            };
        }
    };
}
