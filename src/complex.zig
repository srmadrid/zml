//! Namespace for complex operations.

const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const rational = @import("rational.zig");
const _real = @import("real.zig");

/// Arbitrary-precision complex type, represented as a pair of non-integral
/// allocated numeric values (e.g., `Rational`, `Real`).
pub fn Complex(comptime N: type) type {
    if (!types.isNumeric(N) or !types.isAllocated(N) or types.isIntegral(N))
        @compileError("zml.Complex: N must be a non-integral allocated numeric type, got \n\tN: " ++ @typeName(N) ++ "\n");

    return struct {
        re: N,
        im: N,
        flags: Flags,

        /// Type signature
        pub const is_numeric = true;
        pub const is_complex = true;
        pub const is_complex_type = true;
        pub const is_signed = true;
        pub const is_allocated = true;
        pub const is_custom = types.isCustomType(N);

        /// Scalar type
        pub const Scalar = N;

        pub const empty: Complex(N) = .{
            .re = .empty,
            .im = .empty,
            .flags = .{ .owns_data = false, .writable = false },
        };

        /// Initializes a new complex with the specified real and imaginary
        /// parts.
        ///
        /// ## Signature
        /// ```zig
        /// Complex(N).initSet(allocator: std.mem.Allocator, real: R, imaginary: I) !Complex
        /// ```
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `real` (`anytype`): The value to set the real part to. Must be a
        ///   numeric or a string.
        /// * `imaginary` (`anytype`): The value to set the imaginary part to.
        ///   Must be a numeric or a string.
        ///
        /// ## Returns
        /// `Complex(N)`: The newly initialized complex.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn initSet(allocator: std.mem.Allocator, real: anytype, imaginary: anytype) !Complex(N) {
            var c: Complex(N) = undefined;
            c.re = try .init(allocator, 0, 0);
            c.im = try .init(allocator, 0, 0);
            c.flags = .{ .owns_data = true, .writable = true };
            errdefer c.deinit(allocator);

            try c.set(allocator, real, imaginary);

            return c;
        }

        /// Deinitializes the complex, freeing any allocated memory and
        /// invalidating it.
        ///
        /// If the complex does not own its data, no memory is freed and this
        /// only invalidates it.
        ///
        /// ## Arguments
        /// * `self` (`*Complex(N)`): A pointer to the complex to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *Complex(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                self.re.deinit(allocator);
                self.im.deinit(allocator);
            }

            self.* = undefined;
        }

        /// Sets the value of the complex.
        ///
        /// ## Signature
        /// ```zig
        /// Complex(N).set(allocator: std.mem.Allocator, real: R, imaginary: I) !void
        /// ```
        ///
        /// ## Arguments
        /// * `self` (`*Complex(N)`): A pointer to the complex to set.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `real` (`anytype`): The value to set the real part to. Must be a
        ///   numeric or a string. Complex values use their real part.
        /// * `imaginary` (`anytype`): The value to set the imaginary part to.
        ///   Must be a numeric or a string. Complex values use their real part.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `complex.Error.NotWritable`: If the `Rational` is not writable.
        /// * `complex.Error.NotFinite`: If either value is a float or cfloat
        ///   number that is not finite.
        pub fn set(self: *Complex(N), allocator: std.mem.Allocator, real: anytype, imaginary: anytype) !void {
            // TODO: improve greatly

            if (!self.flags.writable)
                return Error.NotWritable;

            try self.re.set(allocator, real);
            try self.im.set(allocator, imaginary);
        }

        /// Creates a copy of the complex.
        ///
        /// ## Arguments
        /// * `self` (`*const Complex`): A pointer to the complex to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `Complex(N)`: The newly created copy of the complex.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn copy(self: *const Complex(N), allocator: std.mem.Allocator) !Complex(N) {
            var re: N = try self.re.copy(allocator);
            errdefer re.deinit(allocator);
            const im: N = try self.im.copy(allocator);

            return .{
                .re = re,
                .im = im,
                .flags = .{ .owns_data = true, .writable = true },
            };
        }
    };
}

// Arithmetic operations
pub const Add = @import("complex/add.zig").Add;
pub const add = @import("complex/add.zig").add;
pub const add_ = @import("complex/add_.zig").add_;
pub const Sub = @import("complex/sub.zig").Sub;
pub const sub = @import("complex/sub.zig").sub;
pub const sub_ = @import("complex/sub_.zig").sub_;
pub const Mul = @import("complex/mul.zig").Mul;
pub const mul = @import("complex/mul.zig").mul;
pub const mul_ = @import("complex/mul_.zig").mul_;
pub const Div = @import("complex/div.zig").Div;
pub const div = @import("complex/div.zig").div;
pub const div_ = @import("complex/div_.zig").div_;

// Comparison operations
// pub const eq = @import("complex/eq.zig").eq;
// pub const ne = @import("complex/ne.zig").ne;

// Basic operations
// pub const abs = @import("complex/abs.zig").abs;
pub const neg = @import("complex/neg.zig").neg;
pub const conj = @import("complex/conj.zig").conj;

pub const Error = error{
    ZeroDivision,
    NotWritable,
    NotFinite,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    writable: bool = true,
};
