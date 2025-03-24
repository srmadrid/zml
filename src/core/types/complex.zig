const std = @import("std");
const RationalUnmanaged = @import("rational.zig").RationalUnmanaged;

pub const ComplexUnmanaged = struct {
    re: RationalUnmanaged,
    im: RationalUnmanaged,
};

pub const Complex = struct {
    re: RationalUnmanaged,
    im: RationalUnmanaged,
    allocator: std.mem.Allocator,
};
