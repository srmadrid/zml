const std = @import("std");

pub const IntegerUnmanaged = struct {
    limbs: []usize,
    metadata: usize,
};

pub const Integer = struct {
    limbs: []usize,
    metadata: usize,
    allocator: std.mem.Allocator,
};
