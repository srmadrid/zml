const std = @import("std");
const IntegerUnmanaged = @import("integer.zig").IntegerUnmanaged;

pub const RationalUnmanaged = struct {
    numerator: IntegerUnmanaged,
    denominator: IntegerUnmanaged,
};

pub const Rational = struct {
    numerator: IntegerUnmanaged,
    denominator: IntegerUnmanaged,
    allocator: std.mem.Allocator,
};
