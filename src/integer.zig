const std = @import("std");

pub const Integer = struct {
    limbs: []usize,
    size: usize,

    pub const empty: Integer = .{
        .limbs = &.{},
        .size = 0,
    };
};
