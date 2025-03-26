const std = @import("std");
const IntegerUnmanaged = @import("integer.zig").IntegerUnmanaged;
const Integer = @import("integer.zig").Integer;
const RationalUnmanaged = @import("rational.zig").RationalUnmanaged;
const Rational = @import("rational.zig").Rational;
const RealUnmanaged = @import("real.zig").RealUnmanaged;
const Real = @import("real.zig").Real;

// Limit Complex to Rational and Real? i.e., remove Integer? If so, edit coercion functions

pub fn ComplexUnmanaged(comptime T: type) type {
    if (T != IntegerUnmanaged and T != Integer and T != RationalUnmanaged and T != Rational and T != RealUnmanaged and T != Real) @compileError("Unsupported type for ComplexUnmanaged: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,

        pub const empty: ComplexUnmanaged(T) = .{
            .re = .empty,
            .im = .empty,
        };
    };
}

pub fn Complex(comptime T: type) type {
    if (T != IntegerUnmanaged and T != Integer and T != RationalUnmanaged and T != Rational and T != RealUnmanaged and T != Real) @compileError("Unsupported type for Complex: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,
        allocator: std.mem.Allocator,

        pub const empty: Complex(T) = .{
            .re = .empty,
            .im = .empty,
            .allocator = undefined,
        };
    };
}
