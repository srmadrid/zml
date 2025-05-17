const std = @import("std");
const Integer = @import("integer.zig").Integer;

pub const Rational = struct {
    numerator: Integer,
    denominator: Integer,

    pub const empty: Rational = .{
        .numerator = .empty,
        .denominator = .empty,
    };
};
