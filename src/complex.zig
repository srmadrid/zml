const std = @import("std");
const Integer = @import("integer.zig").Integer;
const Rational = @import("rational.zig").Rational;
const Real = @import("real.zig").Real;

// Limit Complex to Rational and Real? i.e., remove Integer? If so, edit coercion functions

pub fn Complex(comptime T: type) type {
    if (T != Integer and T != Rational and T != Real)
        @compileError("Unsupported type for Complex: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,

        pub const empty: Complex(T) = .{
            .re = .empty,
            .im = .empty,
        };
    };
}
