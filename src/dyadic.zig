//! Namespace for dyadic operations.

const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const rational = @import("rational.zig");
const _real = @import("real.zig");

/// Arbitrary-precision (chosen at compile time)
/// dyadic type, represented as two arbitrary-precision
/// rational (`Rational`) or real (`Real`) numbers.
pub fn Dyadic(mantissa_bits: u16, exponent_bits: u16) type {
    if (mantissa_bits == 0 or exponent_bits == 0)
        @compileError("Dyadic requires both mantissa_bits and exponent_bits to be non-zero");

    return packed struct {
        mantissa: @Type(.{
            .int = .{
                .signedness = .unsigned,
                .bits = mantissa_bits,
            },
        }),
        exponent: @Type(.{
            .int = .{
                .signedness = .signed,
                .bits = exponent_bits,
            },
        }),
        positive: bool,

        /// Initializes a new `Dyadic` with the specified value.
        ///
        /// Signature
        /// ---------
        /// ```zig
        /// fn initSet(value: V) Dyadic(mantissa_bits, exponent_bits)
        /// ```
        ///
        /// Parameters
        /// ----------
        /// `value` (`anytype`):
        /// The value to set the dyadic to. Must be a numeric type or a string.
        ///
        /// Returns
        /// -------
        /// `Dyadic`:
        /// The newly initialized `Dyadic`.
        pub fn initSet(value: anytype) Dyadic(mantissa_bits, exponent_bits) {
            _ = value;
            return .{ .mantissa = 0, .exponent = 0, .positive = true }; // TODO: implement
        }

        /// Sets the value of the `Dyadic`.
        ///
        /// Signature
        /// ---------
        /// ```zig
        /// fn set(value: V) void
        /// ```
        ///
        /// Parameters
        /// ----------
        /// `self` (`*Dyadic`):
        /// A pointer to the `Dyadic` to set.
        ///
        /// `value` (`anytype`):
        /// The value to set the dyadic to. Must be a numeric type or a string.
        ///
        /// Returns
        /// -------
        /// `void`
        pub fn set(self: *Dyadic(mantissa_bits, exponent_bits), value: anytype) void {
            // TODO: implement
            _ = value;
            self.* = .{ .mantissa = 0, .exponent = 0, .positive = true };
        }
    };
}

// Arithmetic operations
// pub const add = @import("dyadic/add.zig").add;
// pub const add_ = @import("dyadic/add_.zig").add_;
// pub const sub = @import("dyadic/sub.zig").sub;
// pub const sub_ = @import("dyadic/sub_.zig").sub_;
// pub const mul = @import("dyadic/mul.zig").mul;
// pub const mul_ = @import("dyadic/mul_.zig").mul_;
// pub const div = @import("dyadic/div.zig").div;
// pub const div_ = @import("dyadic/div_.zig").div_;

// Comparison operations
// pub const eq = @import("dyadic/eq.zig").eq;
// pub const ne = @import("dyadic/ne.zig").ne;

// Basic operations
// pub const abs = @import("dyadic/abs.zig").abs;
// pub const neg = @import("dyadic/neg.zig").neg;
