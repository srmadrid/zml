const std = @import("std");
const Integer = @import("integer.zig").Integer;
const IntegerManaged = @import("integer.zig").IntegerManaged;
const Rational = @import("rational.zig").Rational;
const RationalManaged = @import("rational.zig").RationalManaged;
const Real = @import("real.zig").Real;
const RealManaged = @import("real.zig").RealManaged;

// Limit Complex to Rational and Real? i.e., remove Integer? If so, edit coercion functions

pub fn Complex(comptime T: type) type {
    if (T != Integer and T != Rational and T != Real)
        @compileError("Unsupported type for ComplexUnmanaged: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,

        pub const empty: Complex(T) = .{
            .re = .empty,
            .im = .empty,
        };
    };
}

pub fn ComplexManaged(comptime T: type) type {
    if (T != Integer and T != Rational and T != Real)
        @compileError("Unsupported type for ComplexUnmanaged: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,
        allocator: std.mem.Allocator,

        pub const empty: ComplexManaged(T) = .{
            .re = .empty,
            .im = .empty,
            .allocator = undefined,
        };
    };
}
