const std = @import("std");
const IntegerUnmanaged = @import("integer.zig").IntegerUnmanaged;

pub const RationalUnmanaged = struct {
    numerator: IntegerUnmanaged,
    denominator: IntegerUnmanaged,

    pub const empty: RationalUnmanaged = .{
        .numerator = .empty,
        .denominator = .empty,
    };
};

pub const Rational = struct {
    numerator: IntegerUnmanaged,
    denominator: IntegerUnmanaged,
    allocator: std.mem.Allocator,

    pub const empty: Rational = .{
        .numerator = IntegerUnmanaged.empty,
        .denominator = IntegerUnmanaged.empty,
        .allocator = undefined,
    };
};
