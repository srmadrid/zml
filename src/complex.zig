const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const rational = @import("rational.zig");
const _real = @import("real.zig");

/// A complex number type represented as two `Rational` or `Real` numbers.
pub fn Complex(comptime T: type) type {
    if (T != rational.Rational and T != _real.Real)
        @compileError("Unsupported type for Complex: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,
        flags: Flags,

        pub const empty: Complex(T) = .{
            .re = .empty,
            .im = .empty,
            .flags = .{ .owns_data = false, .writable = false },
        };

        /// Initializes a new `Complex` with the specified numerator and
        /// denominator.
        ///
        /// Signature
        /// ---------
        /// ```zig
        /// fn initSet(allocator: std.mem.Allocator, real: R, imaginary: I) !Complex
        /// ```
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `real` (`anytype`):
        /// The value to set the real part to. Must be a numeric type or a
        /// string.
        ///
        /// `imaginary` (`anytype`):
        /// The value to set the imaginary part to. Must be a numeric type or a
        /// string.
        ///
        /// Returns
        /// -------
        /// `Complex`:
        /// The newly initialized `Complex`.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        pub fn initSet(allocator: std.mem.Allocator, real: anytype, imaginary: anytype) !Complex(T) {
            var c: Complex(T) = undefined;
            c.re = try .init(allocator, 0, 0);
            c.im = try .init(allocator, 0, 0);
            c.flags = .{ .owns_data = true, .writable = true };
            errdefer c.deinit(allocator);

            try c.set(allocator, real, imaginary);
            return c;
        }

        /// Deinitializes the `Complex`, freeing any allocated memory and
        /// invalidating it.
        ///
        /// If the `Complex` does not own its data, no memory is freed and this
        /// only invalidates it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*Complex`):
        /// A pointer to the `Complex` to deinitialize.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory deallocation. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// Returns
        /// -------
        /// `void`
        pub fn deinit(self: *Complex(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                self.re.deinit(allocator);
                self.im.deinit(allocator);
            }

            self.* = undefined;
        }

        /// Sets the value of the `Complex`.
        ///
        /// Signature
        /// ---------
        /// ```zig
        /// fn set(allocator: std.mem.Allocator, real: R, imaginary: I) !void
        /// ```
        ///
        /// Parameters
        /// ----------
        /// `self` (`*Complex`):
        /// A pointer to the `Complex` to set.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `real` (`anytype`):
        /// The value to set the real part to. Must be a numeric type or a
        /// string. Complex values use their real part.
        ///
        /// `imaginary` (`anytype`):
        /// The value to set the imaginary part to. Must be a numeric type or a
        /// string. Complex values use their real part.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `complex.Error.NotWritable`:
        /// If the `Rational` is not writable.
        ///
        /// `complex.Error.NotFinite`:
        /// If either value is a float or complex number that is not finite.
        pub fn set(self: *Complex(T), allocator: std.mem.Allocator, real: anytype, imaginary: anytype) !void {
            comptime var S: type = @TypeOf(self);

            comptime if (!types.isPointer(S) or types.isConstPointer(S))
                @compileError("complex.set requires self to be a mutable pointer, got " ++ @typeName(S));

            S = types.Child(S);

            comptime if (types.numericType(S) != .complex)
                @compileError("complex.set requires self to be a complex, got " ++ @typeName(S));

            if (!self.flags.writable)
                return Error.NotWritable;

            try self.re.set(allocator, real);
            try self.im.set(allocator, imaginary);
        }

        /// Creates a copy of the `Complex`.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const Complex`):
        /// A pointer to the `Complex` to copy.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// Returns
        /// -------
        /// `Complex`:
        /// The newly created copy of the `Complex`.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        pub fn copy(self: *const Complex(T), allocator: std.mem.Allocator) !Complex(T) {
            var re: T = try self.re.copy(allocator);
            errdefer re.deinit(allocator);
            const im: T = try self.im.copy(allocator);

            return .{
                .re = re,
                .im = im,
                .flags = .{ .owns_data = true, .writable = true },
            };
        }
    };
}

// Arithmetic operations
pub const add = @import("complex/add.zig").add;
pub const add_ = @import("complex/add_.zig").add_;
pub const sub = @import("complex/sub.zig").sub;
pub const sub_ = @import("complex/sub_.zig").sub_;
pub const mul = @import("complex/mul.zig").mul;
pub const mul_ = @import("complex/mul_.zig").mul_;
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
